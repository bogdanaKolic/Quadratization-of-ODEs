# -*- coding: utf-8 -*-

import random
import math
import collections
import itertools

max_coeff = 100 # maximal generated coefficient of a monomial
max_power = 100 # maximal generated power of a variable in a monomial
max_width = 100 # maximal number of original variables
max_length = 50 # maximal number of monomials a generated equation can have

# dictionary relating certain systems with their optimal monomial quadratization
optimal_quadratization = {
    'monom3' : {(0, 2, 1), (1, 2, 0), (3, 0, 0), (0, 3, 0), (0, 0, 3), (1, 2, 2),
             (2, 1, 2), (1, 0, 2), (1, 1, 1), (2, 2, 1)},
    'circular(2,3)' : {(1, 1), (2, 0), (0, 2)},
    'cubcic_cycle(3)' : {(3, 0, 0), (0, 3, 0), (0, 0, 3), (2, 0, 0), (0, 2, 0), (0, 0, 2)},
    'cubic_cycle(7)' : {(3, 0, 0, 0, 0, 0, 0), (0, 3, 0, 0, 0, 0, 0), (0, 0, 3, 0, 0, 0 ,0),
             (0, 0, 0, 3, 0, 0, 0), (0, 0, 0, 0, 3 ,0 ,0), (0, 0, 0, 0, 0, 3 ,0),
             (0, 0, 0, 0, 0, 0, 3), (2, 0, 0, 0, 0 ,0 ,0), (0, 2, 0, 0, 0, 0 ,0), 
             (0, 0, 2, 0, 0 ,0 ,0), (0, 0, 0, 2, 0, 0, 0), (0, 0, 0, 0, 2, 0, 0),
             (0, 0, 0, 0, 0, 2 ,0), (0 ,0, 0, 0, 0, 0, 2)}
    }

class Monomial:
    """ Monomials are represented with a coefficient, number of variables, 
    a tuple of powers for each variable, boolean variables to indicate their
    presence in the system and (if it exists) their origin - the substitution
    whose derivative they belong to"""
    def __init__(self, n = 0, c = 1, var = tuple(), origin = None):
        self.width = n # number of variables
        self.coefficient = c
        self.variables = var
        self.is_currently_present = False # does this monomial appear in the 
                                          # system with currently considered 
                                          # substitutions
        self.is_always_present = False # is this monomial a part of the
                                       # original system
        self.origin = origin # if origin is None, the monomial is not a part of 
                             # a derivative or it belongs to the original system
                             # of ODEs, otherwise the origin is equal to the
                             # substituion whose derivative this monomial 
                             # belongs to
        
    def set_random_data(self, n = 0):
        self.width = n # number of variables
        self.coefficient = random.randint(1, max_coeff)
        self.variables = tuple()
        for i in range(n):
            self.variables += (random.randint(0, max_power), )
        
    def copy(self):
        mono = Monomial()
        mono.coefficient = self.coefficient
        mono.width = self.width
        mono.variables = self.variables
        return mono
    
    def change_degree(self, index):
        """ Auxiliary function for calculating a derivative - it takes a monomial 
        and reduces the power of variable 'index' by 1 and calculates the
        new coefficient"""
        coeff = self.coefficient
        t = tuple()
        for i in range(self.width):
            if i == index:
                coeff *= self.variables[i]
                t += (self.variables[i] - 1, )
            else:
                t += (self.variables[i], )
        self.variables = t
        self.coefficient = coeff
    
    def calculate_derivative(self, equations):
        """ Calculates a polynomial that represents the derivative of a monomial
        given the list of derivatives of the original variables xi (equations)"""
        poly = Polynomial()
        for i, deg in enumerate(self.variables):
            if deg != 0:
                mono = self.copy()
                mono.change_degree(i)
                poly.add_polynomial(Polynomial.multiply(equations[i], mono))
        poly.shrink()
        return poly
    
    def multiply(monom_1, monom_2):
        """ Multiplies two given monomials and return a resulting polynomial"""
        coeff = monom_1.coefficient * monom_2.coefficient
        product_variables = tuple([x + y for (x, y) in zip(monom_1.variables, monom_2.variables)])
        product = Monomial(monom_1.width, coeff, product_variables)
        return product
    
    def __add__(self, other):
        if not (self.origin is None):
            self.is_currently_present = self.origin.is_currently_present
        if not (other.origin is None):
            other.is_currently_present = other.origin.is_currently_present
        return self.is_currently_present or other.is_currently_present
    
    def __mul__(self, other):
        if not (self.origin is None):
            self.is_currently_present = self.origin.is_currently_present
        if not (other.origin is None):
            other.is_currently_present = other.origin.is_currently_present
        return self.is_currently_present and other.is_currently_present
    
    def __eq__(self, other):
        return self.variables == other.variables
    
    def __hash__(self):
        return hash(self.variables)
    
    def __repr__(self):
        return str(self.coefficient) + '*' + str(self.variables)      
    
    def __str__(self, use_coefficient = True):
        if use_coefficient:
            s = str(self.coefficient)
        else:
            s = ''
        if self.variables != (0, ) * self.width:
            s += ' * ('
            for i, power in enumerate(self.variables):
                if power > 0:
                    s += f'x{i}^{power} * '
                elif power < 0:
                    s += f'x{i}^({power}) * '
            s = s.strip(' *')
            s += ')'
        if s == '':
            s = '1'
        return s
    
class Substitution(Monomial):
    """Substituion monomials in addition have a pre-calculated derivative and a 
    number of quadratic monomials they are factors of; they have allowed 
    negative powers(Laurent monomials)"""
    def __init__(self, monomial, equations, index = -1):
        self.number_of_quadratic_monomials = 0 # the number of different quadratic  
                                               # monomials it is a factor of
        if index == -1:
            Monomial.__init__(self, len(monomial), 1, monomial)
        else:
            Monomial.__init__(self, monomial.width, monomial.coefficient, tuple())
            for i, deg in enumerate(monomial.variables):
                if i == index:
                    self.variables += (monomial.variables[i] - 1, )
                else:
                    self.variables += (monomial.variables[i], )
        self.derivative = self.calculate_derivative(equations)
        for monom_derivative in self.derivative.monomials:
            monom_derivative.origin = self
        #self.index = index
    
    def get_number_of_quadratic_monomials(substitution):
        return substitution.number_of_quadratic_monomials
    
    def multiply_x(self, i):
        """Multiplies Laurent monomial by xi"""
        product  = Monomial(self.width, self.coefficient)
        l = [x for x in self.variables]
        l[i] += 1
        product.variables = tuple(l)
        return product
    
    def __str__(self):
        s = super().__str__(False)
        #s += '\nderivative: ' + str(self.derivative) + '\nindex: ' + str(self.index)
        return s

class Polynomial:
    
    def __init__(self, n = 0, m = 0, monomials =[]):
        self.width = n #number of variables in a monomial
        self.length = m # number of monomials in a polynomial
        self.monomials = monomials.copy()
        if len(self.monomials) != m:
            for i in range(m - len(self.monomials)):
                mono = Monomial(n)
                self.monomials.append(mono)
    
    def set_random_data(self, n = 0, m = 0):
        self.width = n #number of variables in a monomial
        self.length = m # number of monomials in a polynomial
        self.monomials = []
        for i in range(m):
            mono = Monomial()
            mono.set_random_data(n)
            self.monomials.append(mono)
    
    def add_monomial(self, mono):
        """ Adds a monomial to a polynomial """
        self.monomials.append(mono)
        self.length += 1
        if self.width == 0:
            self.width = mono.width
    
    def add_polynomial(self, poly):
        """ Adds a polynomial to a polynomial """
        for mono in poly.monomials:
            self.add_monomial(mono)
    
    def multiply(poly, monom):
        """ Multiplies a polynomial by a monomial and return the product """
        product = Polynomial()
        for mono in poly.monomials:
            product.add_monomial(mono.multiply(monom))
        return product
    
    def shrink(self):
        """ Compress the polynomial into a simpler form """
        monomials_dictionary = {}
        for monom in self.monomials:
            if monom in monomials_dictionary:
                monomials_dictionary[monom] += monom.coefficient
            else:
                monomials_dictionary[monom] = monom.coefficient
        monomials = []
        for mono in monomials_dictionary:
            mono.coefficient = monomials_dictionary[mono]
            if mono.coefficient != 0:
                monomials.append(mono)
        self.monomials = monomials.copy()
    
    def __str__(self):
        s = ''
        for i in range(self.length - 1):
            s += str(self.monomials[i]) + ' + '
        if self.monomials != []:
            s += str(self.monomials[- 1])
        return s
    
class Equation(Polynomial):
    """ A class for representing given differential equations (derivatives of
    the original x variables) """
    def __init__(self, index, n = 0, m = 0, monomials = []):
        Polynomial.__init__(self, n, m, monomials)
        self.index = index # corresponds to the variable being differentiated
        
    def set_random(self, n = 0):
        self.width = n
        m = random.randint(1, max_length)
        self.length = m
        self.set_random_data(n, m)
        
    def calculate_substitutions(self, equations):
        """ Calculates all the Laurent monomials that are obtained from one 
        equation """
        substitutions = set()
        for m in self.monomials:
            substitutions.add(Substitution(m, equations, self.index))
        return substitutions
    
    def __str__(self):
        s = f'x{self.index} = ' + super().__str__()
        return s
    
class QuadratizationProblem():
    """ Class for performing a single test - either on a benchmark set or 
    generates a random set of equations and tries to minimize the number of 
    Laurent monomials used in quadratization """
    def __init__(self, width = 0, equations = [], additional_sub = set()):
        self.width = width # number of original variables
        self.equations = equations.copy() # the original system
        self.all_substitutions = set() # set of all substitutions considered 
                                       # for addition
        self.optional_substitutions = [] # list of substitutions that may or 
                                         # may not be present in the quadratized system
        self.original_substitutions = set() # set of Laurent monomials from 
                                            # the proposition
        self.optimal_solution = set()
        self.current_length = 0 # number of substitutions currently being 
                                # considered for quadratization
        self.min_length = math.inf # minimal length of a quadratization
        self.product_results = {} # the dictionary connecting the products to
                                  # their factors
        self.at_most_quadratic_monomials = {} # the inverse dictionary connecting
                                              # substitutions to the monomials 
                                              # they are factors of
                                              
        self.additional_substitutions = additional_sub.copy() # set of substitutions
                                                              # added to the orginal set 
                                                              # mentioned in the propositon
        self.precompute_part1()
    
    def precompute_part1(self):
        """ Precomputation before the additional variables are computed"""
        # compute the set of Laurent monomials from the proposition
        for e in self.equations:
            self.original_substitutions = self.original_substitutions.union(e.calculate_substitutions(self.equations))
        # remove unnecessary substitutions
        set_of_original = set()
        for o in self.original_substitutions:
            powers = sum([abs(x) for x in o.variables])
            if powers > 1:
                set_of_original.add(o)
            if sum(o.variables) == -1:
                set_of_original.add(o)
        self.original_substitutions = set_of_original.copy()
        self.min_length = len(self.original_substitutions)
        # add the known solution
        self.optimal_solution = {y for y in self.original_substitutions}
    
    def precompute_part2(self):
        """ Precomputation after the additional variables are computed"""
        # remove unnecessary substitutions
        set_of_additional = set()
        for a in self.additional_substitutions:
            powers = sum([abs(x) for x in a.variables])
            if powers > 1:
                set_of_additional.add(a)
            if sum(a.variables) == -1:
                    set_of_additional.add(a)
        self.additional_substitutions = set_of_additional.copy()
        self.all_substitutions = self.original_substitutions.union(self.additional_substitutions)
        # pre-compute the dictionary of products represented as at most quadratic
        # monomials appearing in the system with all the substitutions
        search_space_deg1 = set()
        # add 'a constant'
        constant_monom = Monomial(self.width, 1, (0, ) * self.width)
        constant_monom.is_currently_present = True
        constant_monom.is_always_present = True
        search_space_deg1.add(constant_monom)
        # add variables xi
        l = [0] * self.width
        for i in range(self.width):
            l[i] = 1
            monom_x = Monomial(self.width, 1, tuple(l))
            monom_x.is_currently_present = True
            monom_x.is_always_present = True
            search_space_deg1.add(monom_x)
            l[i] = 0
        # add substitution variables
        search_space_deg1 = search_space_deg1.union(self.all_substitutions)
        search_space_deg1 = list(search_space_deg1)
        # compute the set of monomials appearing in the derivatives
        monomials_in_derivatives = {}
        for e in self.equations:
            for monom in e.monomials:
                if monom.variables not in monomials_in_derivatives:
                    monomials_in_derivatives[monom.variables] = [monom]
                else:
                    monomials_in_derivatives[monom.variables].append(monom)
        for substitution in self.all_substitutions:
            for monom_derivative in substitution.derivative.monomials:
                if monom_derivative.variables not in monomials_in_derivatives:
                    monomials_in_derivatives[monom_derivative.variables] = [monom_derivative]
                else:
                    monomials_in_derivatives[monom_derivative.variables].append(monom_derivative)
        # compute the at most quadratic monomials
        for i in range(len(search_space_deg1)):
            for j in range(i, len(search_space_deg1)):
                first_monom = search_space_deg1[i]
                second_monom = search_space_deg1[j]
                product = Monomial.multiply(first_monom, second_monom)
                if product.variables in monomials_in_derivatives:
                    for p in monomials_in_derivatives[product.variables]:
                        if first_monom not in self.at_most_quadratic_monomials:
                            self.at_most_quadratic_monomials[first_monom] = {p}
                        else:
                            self.at_most_quadratic_monomials[first_monom].add(p)
                        if second_monom not in self.at_most_quadratic_monomials:
                            self.at_most_quadratic_monomials[second_monom] = {p}
                        else:
                            self.at_most_quadratic_monomials[second_monom].add(p)
                        if p not in self.product_results:
                            self.product_results[p] = {(first_monom, second_monom)}
                        else:
                            self.product_results[p].add((first_monom, second_monom))
            
    def find_must_have_substitutions(self):
        """ Find out which substitutions have to be present in the system 
        and compute the list of optional substitutions"""
        check = collections.deque()
        for monom in self.product_results:
            if monom.origin is None:
                check.append(monom)
        while len(check) != 0:
            monom = check.popleft()
            if len(self.product_results[monom]) == 1:
                for element in self.product_results[monom]:
                    for m in element:
                        if not m.is_always_present and m.origin is not None:
                            check.append(m)
                        m.is_currently_present = True
                        m.is_always_present = True
        for substitution in self.all_substitutions:
            if substitution.is_always_present:
                self.current_length += 1
            else:
                self.optional_substitutions.append(substitution)
        for substitution in self.all_substitutions:
            if substitution in self.at_most_quadratic_monomials:
                self.at_most_quadratic_monomials[substitution] = {product for product in self.at_most_quadratic_monomials[substitution] if product.origin is None or product.origin in self.all_substitutions}
                substitution.number_of_quadratic_monomials = len(self.at_most_quadratic_monomials[substitution])
        self.optional_substitutions = sorted(self.optional_substitutions, key = Substitution.get_number_of_quadratic_monomials)

    def remove_impossible_substitutions(self):
        """ Find out which substitutions cannot be present in the system"""
        while(True):
            to_remove = set()
            for substitution in self.additional_substitutions:
                if substitution not in self.at_most_quadratic_monomials or len(self.at_most_quadratic_monomials[substitution]) == 0:
                    to_remove.add(substitution)
                    continue
                for monom in substitution.derivative.monomials:
                    if (monom not in self.product_results) or len(self.product_results[monom]) == 0:
                        to_remove.add(substitution)
            for substitution in to_remove:
                if substitution not in self.original_substitutions:
                    if substitution in self.at_most_quadratic_monomials:
                        for product in self.at_most_quadratic_monomials[substitution]:
                            self.product_results[product] = {factors for factors in self.product_results[product] if (factors[0] != substitution and factors[1] != substitution)}
                    self.additional_substitutions.remove(substitution)
                    self.all_substitutions.remove(substitution)
            if len(to_remove) == 0:
                return
        
    def random_test(self):
        self.width = random.randint(1, max_width)
        print('width: ', self.width)
        self.equations = []
        for i in range(self.width):
            eq = Equation(i, self.width)
            eq.set_random(self.width)
            self.equations.append(eq)
        self.precompute_part1()
        self.precompute_part2()
        self.remove_impossible_substitutions()
        self.find_must_have_substitutions()
    
    def load_from_file(self, filename):
        with open(filename) as infile:
            index  = 0
            self.width = int(infile.readline())
            for line in infile:
                monomials = line.split('+')
                eq = []
                for m in monomials:
                    index_of_star = m.find('*')
                    if index_of_star != -1:
                        coef = int(m[: index_of_star])
                        var = m[index_of_star + 1: ]
                    else:
                        coef = 1
                        var = m
                    variables = tuple()
                    var = var.strip('( )\n')
                    var = var.split(',')
                    for power in var:
                        variables += (int(power), )
                    eq.append(Monomial(self.width, coef, variables))
                self.equations.append(Equation(index, self.width, len(eq), eq))
                index += 1
        self.precompute_part1()
        
    def precompute_non_quadratized(self):
        """ Retrieve the set of non-quadratized monomials present in the system """
        non_quadratized = set()
        # add the monomials from the original equations that aren't quadratized
        for eq in self.equations:
            for monom in eq.monomials:
                result = False
                for factors in self.product_results[monom]:
                    result = result or factors[0] * factors[1]
                    if result:
                        break
                if not result:
                    non_quadratized.add(monom)
        # add the monomials from the must-have substitutions
        for substitution in self.all_substitutions:
            if substitution.is_currently_present:
                for monom_derivative in substitution.derivative.monomials:
                    result = False
                    for factors in self.product_results[monom_derivative]:
                        result = result or factors[0] * factors[1]
                        if result:
                            break
                    if not result:
                        non_quadratized.add(monom_derivative)
        return non_quadratized
    
    def quadratizes(self, substitution, non_quadratized):
        """Given a substitution and a set of non_quadratized monomials in the
        system without the substitution, computes a new set of non_quadratized
        monomials in a system with the substitution"""
        still_not_quadratized = set()
        for monom in non_quadratized:
            result = False
            for factors in self.product_results[monom]:
                result = result or factors[0] * factors[1]
                if result:
                    break
            if not result:
                still_not_quadratized.add(monom)
        for monom_derivative in substitution.derivative.monomials:
            result = False
            for factors in self.product_results[monom_derivative]:
                result = result or factors[0] * factors[1]
                if result:
                    break
            if not result:
                still_not_quadratized.add(monom_derivative)
        return still_not_quadratized
    
    def variables_from_derivatives(self):
        """ Calculate the additional variables by dividing each monomial 
        appearing in the derivatives of the original substitutions by 
        original variables """
        for substitution in self.original_substitutions:
             for monom_derivative in substitution.derivative.monomials:
                 for i in range(self.width):
                     self.additional_substitutions.add(Substitution(monom_derivative, self.equations, i))
    
    def random_new_variables_from_derivatives(self):
        """ Randomly compute substitutions based on the derivatives of the 
        original substitutions"""
        for substitution in self.original_substitutions:
             for monom_derivative in substitution.derivative.monomials:
                 l = [random.randint(-1, x + 1) for x in monom_derivative.variables]
                 p = [x - y for (x, y) in zip(monom_derivative.variables, l)]
                 self.additional_substitutions.add(Substitution(tuple(l), self.equations))
                 self.additional_substitutions.add(Substitution(tuple(p), self.equations))
        for equation in self.equations:
             for monom in equation.monomials:
                 l = [random.randint(-1, x) for x in monom.variables]
                 p = [x - y for (x, y) in zip(monom.variables, l)]
                 self.additional_substitutions.add(Substitution(tuple(l), self.equations))
                 self.additional_substitutions.add(Substitution(tuple(p), self.equations))
    
    def variables_derivatives_roots(self):
        """ Additional substitutions calculated as the square roots of the 
        monomials appearing in the derivatives of the original substitutions"""
        for substitution in self.original_substitutions:
             for monom_derivative in substitution.derivative.monomials:
                 l = [x // 2 for x in monom_derivative.variables]
                 self.additional_substitutions.add(Substitution(tuple(l), self.equations))
                 l = [(x + 1) // 2 for x in monom_derivative.variables]
                 self.additional_substitutions.add(Substitution(tuple(l), self.equations))   
      
    def variables_constant_degree(self, degree):
        """ Function for creating additional substitution variables of a 
        constant degree"""
        l = [i for i in itertools.combinations_with_replacement(range(self.width), degree)]
        for comb in l:
            substitution_variables = [0] * self.width
            for x in comb:
                substitution_variables[x] += 1
            self.additional_substitutions.add(Substitution(tuple(substitution_variables), self.equations))
        """l = [j for j in itertools.combinations_with_replacement([i for i in range(-1, 4)], n)]
        for comb in l: 
            for s in itertools.permutations(comb):
                substitutions.add(quadratization.Substitution(s, equations))
                #substitutions.add(s)"""
    
    def reduce(self, position, non_quadratized):
        """ Find out how many Laurent monomials can be neglected by using 
        a recursion on a set of optional substitutions and considering how the 
        substitution at a given position changes the system"""
        if len(non_quadratized) == 0:
            if self.current_length < self.min_length:
                self.min_length = self.current_length
                self.optimal_solution = {y for y in self.optional_substitutions if y.is_currently_present}
        elif self.current_length >= self.min_length or position == len(self.optional_substitutions):
            return None
        else:
            self.optional_substitutions[position].is_currently_present = False
            self.reduce(position + 1, non_quadratized)
            self.optional_substitutions[position].is_currently_present = True
            self.current_length += 1
            new_set = self.quadratizes(self.optional_substitutions[position], non_quadratized)
            self.reduce(position + 1, new_set)
            self.optional_substitutions[position].is_currently_present = False
            self.current_length -= 1
    
    def run(self, strategy, degree_sub = 0, optquad = None):
        """ Search for the quadratization """
        # print(self)
        if strategy == 'from derivatives':
            self.variables_from_derivatives()
        elif strategy == 'random from derivatives':
            self.random_new_variables_from_derivatives()
        elif strategy == 'constant degree':
            self.variables_constant_degree(degree_sub)
        elif strategy == 'roots of monom derivatives':
            self.variables_derivatives_roots()
        elif strategy == 'optimal quadratization':
            for s in optimal_quadratization[optquad]:
                self.additional_substitutions.add(Substitution(s, self.equations))
        self.precompute_part2()
        self.remove_impossible_substitutions()
        self.find_must_have_substitutions()
        non_quadratized = self.precompute_non_quadratized() # non-quadratized monomials
        # look for the solution
        self.reduce(0, non_quadratized)
        for substitution in self.all_substitutions:
            if substitution.is_always_present:
                self.optimal_solution.add(substitution)
        # print('Solution:')
        # for laurent in self.optimal_solution:
        #    print(laurent)
        print(strategy)
        print('number of must-have substitutions: ', len(self.all_substitutions) - len(self.optional_substitutions))
        print('optimal number of substitutions: ', self.min_length)
        print('number of original substitutions: ', len(self.original_substitutions))
        print('number of all substitutions: ', len(self.all_substitutions), '\n')
       
    def __str__(self):
        # s = 'width = ' + str(self.width) + '\n'
        # s += 'original differential equations:\n'
        # for e in self.equations:
        #     s += str(e) + '\n'
        s = 'original substitutions:\n'
        for a in self.original_substitutions:
            s += str(a) + '\n'
        return s       

def main_random(strategy = 'no additional variables', degree_sub = 0, optquad = None):
    """ Simulate quadratizing a set of randomly generated equations """
    for _ in range(10):
        t = QuadratizationProblem()
        t.random_test()
        t.run(strategy, degree_sub, optquad)  

def main_from_file(filename, strategy = 'no additional variables', degree_sub = 0, optquad = None):
    """ Perform a test on the system stored in a file"""
    t = QuadratizationProblem()
    t.load_from_file(filename)
    t.run(strategy, degree_sub, optquad)
    # print(t.min_length)
    for laurent in t.optimal_solution:
        print(laurent.variables)
