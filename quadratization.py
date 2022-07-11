# -*- coding: utf-8 -*-

import random
import math
import collections

max_coeff = 100 # maximal generated coefficient of a monomial
max_power = 100 # maximal generated power of a variable in a monomial
max_width = 100 # maximal number of original variables
max_length = 50 # maximal number of monomials a generated equation can have

class Monomial:
    """ Monomials are represented with a coefficient, number of variables, 
    and a tuple of the powers of each variable """
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
                             # Laurent monomial substituion whose derivative 
                             # this monomial belongs to
        
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
        given the list of derivatives of original x variables (equations)"""
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
    
class Laurent(Monomial):
    """Laurent monomials in addition have the index of the variable by which 
    they are divided and a pre-calculated derivative"""
    def __init__(self, monomial, index, equations):
        Monomial.__init__(self, monomial.width, monomial.coefficient, tuple())
        for i, deg in enumerate(monomial.variables):
            if i == index:
                self.variables += (monomial.variables[i] - 1, )
            else:
                self.variables += (monomial.variables[i], )
        self.derivative = self.calculate_derivative(equations)
        for monom_derivative in self.derivative.monomials:
            monom_derivative.origin = self
        self.index = index
    
    
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
    the x variables) """
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
        substitutions = []
        for m in self.monomials:
            substitutions.append(Laurent(m, self.index, equations))
        return substitutions
    
    def __str__(self):
        s = f'x{self.index} = ' + super().__str__()
        return s
    
class QuadratizationProblem():
    """ Class for performing single test - either on a benchmark set or 
    generates a random set of equations and tries to minimize the number of 
    Laurent monomials used in quadratization """
    def __init__(self, width = 0, equations = []):
        self.width = width
        self.equations = equations.copy()
        self.all_substitutions = []
        self.optional_substitutions = []
        self.optimal_solution = []
        self.current_length = 0
        self.min_length = math.inf
        self.product_results = {} # the dictionary connecting the factors to
                                  # their products
        self.precompute()
        self.find_must_have_substitutions()
        
    def precompute(self):
        for e in self.equations:
            self.all_substitutions += e.calculate_substitutions(self.equations)
        self.all_substitutions = list(set(self.all_substitutions))
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
        for substitution in self.all_substitutions:
            search_space_deg1.add(substitution)
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
                        if p not in self.product_results:
                            self.product_results[p] = [(first_monom, second_monom)]
                        else:
                            self.product_results[p].append((first_monom, second_monom))
    
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
                for m in self.product_results[monom][0]:
                    if not m.is_always_present and m.origin is not None:
                        check.append(m)
                    m.is_currently_present = True
                    m.is_always_present = True
        for substitution in self.all_substitutions:
            if substitution.is_always_present:
                self.current_length += 1
            else:
                self.optional_substitutions.append(substitution)
            
    def random_test(self):
        self.width = random.randint(1, max_width)
        print('width: ', self.width)
        self.equations = []
        for i in range(self.width):
            eq = Equation(i, self.width)
            eq.set_random(self.width)
            self.equations.append(eq)
        self.precompute()
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
        self.precompute()
        self.find_must_have_substitutions()
        
    def precompute_non_quadratized(self):
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
    
    def reduce(self, position, non_quadratized):
        """ Find out how many Laurent monomials can be neglected by using 
        a recursion on a set of optional substitutions and considering how the 
        substitution at a given position changes the system"""
        if len(non_quadratized) == 0:
            if self.current_length < self.min_length:
                self.min_length = self.current_length
                self.optimal_solution = [y for y in self.optional_substitutions if y.is_currently_present]
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
    
    def run(self):
        """ Search for the quadratization """
        # print(self)
        non_quadratized = self.precompute_non_quadratized() # non-quadratized monomials
        self.reduce(0, non_quadratized)
        for substitution in self.all_substitutions:
            if substitution.is_always_present:
                self.optimal_solution.append(substitution)
        # print('Solution:')
        # for laurent in self.optimal_solution:
        #    print(laurent)
        print('optimal number of substitutions: ', self.min_length)
        print('number of all substitutions: ', len(self.all_substitutions), '\n')
        
    def __str__(self):
        s = 'width = ' + str(self.width) + '\n'
        s += 'original differential equations:\n'
        for e in self.equations:
            s += str(e) + '\n'
        #s += 'substitutions:\n'
        #for a in self.all_substitutions:
        #    s += str(a) + '\n'
        return s       

def main_random():
    """ Simulate quadratizing a set of randomly generated equations """
    for _ in range(10):
        t = QuadratizationProblem()
        t.random_test()
        t.run()  

def main_from_file():
    """ Perform a test on the system stored in a file"""
    t = QuadratizationProblem()
    t.load_from_file('cubic_bicycle(7).txt')
    t.run()
