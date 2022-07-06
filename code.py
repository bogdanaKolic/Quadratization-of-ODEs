# -*- coding: utf-8 -*-

import random
import math

max_coeff = 100 # maximal generated coefficient of a monomial
max_power = 30 # maximal generated power of a variable in a monomial
max_width = 3 # maximal number of original variables
max_length = 7 # maximal number of monomials a generated equation can have

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
        self.origin = origin # if None, the monomial is a part of the original
                             # system, otherwise it is a part of the Laurent
                             # monomial substitution
        
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
                medjuproizvod = Polynomial.multiply(equations[i], mono)
                poly.add_polynomial(medjuproizvod)
        poly.shrink()
        return poly
    
    def multiply(monom_1, monom_2):
        """ Multiplies two given monomials and return a resulting polynomial"""
        coeff = monom_1.coefficient * monom_2.coefficient
        product_variables = tuple([x + y for (x, y) in zip(monom_1.variables, monom_2.variables)])
        product = Monomial(monom_1.width, coeff, product_variables)
        return product
    
    def __add__(self, other):
        if self.origin != None:
            self.is_currently_present = self.origin.is_currently_present
        if other.origin != None:
            other.is_currently_present = other.origin.is_currently_present
        return self.is_currently_present or other.is_currently_present
    
    def __mul__(self, other):
        if self.origin != None:
            self.is_currently_present = self.origin.is_currently_present
        if other.origin != None:
            other.is_currently_present = other.origin.is_currently_present
        return self.is_currently_present and other.is_currently_present
    
    def __hash__(self):
        return hash(self.variables)
    
    def __str__(self):
        return str(self.coefficient) + '*' + str(self.variables)      
    
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
        s = super().__str__()
        s += '\nderivative: ' + str(self.derivative) + '\nindex: ' + str(self.index)
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
    
class Test():
    """ Class for performing single test - either on a benchmark set or 
    generates a random set of equations and tries to minimize the number of 
    Laurent monomials used in quadratization """
    def __init__(self, width = 0, equations = []):
        self.width = width
        self.equations = equations.copy()
        self.all_substitutions = []
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
        self.optimal_solution = [y for y in self.all_substitutions]
        self.min_length = len(self.optimal_solution)
        # pre-compute the dictionaries of at most quadratic monomials and the
        # products they represent
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
        """ Find out which substitutions have to be present in the system """
        check = {monom for monom in self.product_results if monom.origin == None}
        newly_found_substitutions = set()
        while len(check) != 0:
            for monom in check:
                if len(self.product_results[monom]) == 1:
                    monom_1 = self.product_results[monom][0][0]
                    monom_2 = self.product_results[monom][0][1]
                    if not monom_1.is_always_present and monom_1.origin != None:
                        newly_found_substitutions.add(monom_1)
                    if not monom_2.is_always_present and monom_2.origin != None:
                        newly_found_substitutions.add(monom_2)
                    monom_1.is_currently_present = True
                    monom_1.is_always_present = True
                    monom_2.is_currently_present = True
                    monom_2.is_always_present = True
            check = newly_found_substitutions.copy()
            newly_found_substitutions.clear()
            
    def random_test(self):
        self.width = random.randint(1, max_width)
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
        
    def is_quadratization(self):
        """ Checks if the set of currently present substitutions works as 
        a quadratization """
        # check whether all the original equations can be quadratized
        for eq in self.equations:
            for monom in eq.monomials:
                result = False
                for factors in self.product_results[monom]:
                    result = result or factors[0] * factors[1]
                if not result:
                    return False
        # check whether all the substitutions used can be quadratized
        for substitution in self.all_substitutions:
            for monom_derivative in substitution.derivative.monomials:
                result = not substitution.is_currently_present
                for factors in self.product_results[monom_derivative]:
                    result = result or factors[0] * factors[1]
                if not result:
                    return False
        return True

    def reduce(self, position):
        """ Find out how many Laurent monomials can be neglected by using 
        a recursion on a set of all substitutions and considering the monomial 
        at a given position"""
        if position == len(self.all_substitutions):
            if self.current_length < self.min_length and self.is_quadratization():
                self.min_length = self.current_length
                self.optimal_solution = [y for y in self.all_substitutions if y.is_currently_present]
        elif self.current_length >= self.min_length:
            return None
        else:
            if not self.all_substitutions[position].is_always_present:
                self.all_substitutions[position].is_currently_present = False
                self.reduce(position + 1)
            self.all_substitutions[position].is_currently_present = True
            self.current_length += 1
            self.reduce(position + 1)
            self.current_length -= 1
    
    def run(self):
        """ Search for the quadratization """
        print(self)
        self.reduce(0)
        print('Solution:')
        for laurent in self.optimal_solution:
           print(laurent.variables)
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
    """ Simulate quadratizing a set of equations """
    for _ in range(10):
        t = Test()
        t.random_test()
        t.run()
    
def main_from_file():
    """ Perform a test on the system stored in a file"""
    t = Test()
    t.load_from_file('monom3.txt')
    t.run()

def generate_circular_test():
    """Create an instance of class Test that represents a circular system 
    of ODEs"""
    width = random.randint(1, 15) # n
    degree = random.randint(1, 50) # k
    equations = []
    l = [0]*width
    for i in range(width):
        l[(i + 1) % width] = degree
        monom = Monomial(width, 1, tuple(l))
        equations.append(Equation(i, width, 1, [monom]))
        l[(i + 1) % width] = 0
    return Test(width, equations)
    

def circular_benchmark_tests(repeat):
    """Perform a series of tests on circular systems of ODEs, for a given
    number of repetitions"""
    with open('circular_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_circular_test()
            outfile.write(str(test))
            test.run()
            outfile.write('\nSoultion:\n')
            for laurent in test.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\n\n')
            
def generate_hardk_test():
    """Create an instance of class Test that represents a hardk system 
    of ODEs, return it and the degree k of the system"""
    width = 3 # variables A(index 0), B(index 1) and C(index 2)
    degree = random.randint(1, 50) # k
    monomials_of_A = [Monomial(width, 1, (0, 0, degree)), Monomial(width, 1, (2, 2, degree - 1))]
    eqA = Equation(0, width, 2, monomials_of_A)
    monomials_of_B = [Monomial(width, 1, (2, 0, 0))]
    eqB = Equation(1, width, 1, monomials_of_B)
    monomials_of_C = [Monomial(width, 1, (0, 2, 0))]
    eqC = Equation(2, width, 1, monomials_of_C)
    return Test(width, [eqA, eqB, eqC]), degree

def hardk_benchmark_tests(repeat):
    """Perform a series of tests on hardk systems of ODEs, for a given
    number of repetitions"""
    with open('hardk_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test, k = generate_hardk_test()
            outfile.write(f'Degree k: {k}\n')
            outfile.write(str(test))
            test.run()
            outfile.write('\nSolution:\n')
            for laurent in test.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\n\n')
            
def generate_monomn_test():
    """Create an instance of class Test that represents a monomn system 
    of ODEs"""
    width = random.randint(1, 9)
    equations = []
    common_monomial = Monomial(width, 1, (2, ) * width)
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 2
        monom = Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(Equation(i, width, 2, [common_monomial, monom]))
    return Test(width, equations)

def monomn_benchmark_tests(repeat):
    """Perform a series of tests on monomn systems of ODEs, for a given
    number of repetitions"""
    with open('monom_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_monomn_test()
            outfile.write(str(test))
            test.run()
            outfile.write('\nSolution:\n')
            for laurent in test.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\n\n')

def generate_hilln_tests():
    """Create two instances of class Test that represent hill systems
    of ODEs, one with 3, and one with 4 variables; return them and their order"""
    order = random.randint(1, 100) # n
    width = 3 # variables H(index 0), I(index 1) and T(index 2)
    monom_H_3 = Monomial(width, order, (0, 2, order - 1))
    monom_I_3 = Monomial(width, - order, (0, 2, order - 1))
    monom_T_3 = Monomial(width, 1, (0, 0, 0))
    equation_H_3 = Equation(0, width, 1, [monom_H_3])
    equation_I_3 = Equation(1, width, 1, [monom_I_3])
    equation_T_3 = Equation(2, width, 1, [monom_T_3])
    test_3 = Test(width, [equation_H_3, equation_I_3, equation_T_3])
    width = 4 # variables H(index 0), I(index 1) and T(index 2) and X(index 4)
    monom_H_4 = Monomial(width, order, (0, 2, order - 1, 1))
    monom_I_4 = Monomial(width, - order, (0, 2, order - 1, 1))
    monom_T_4 = Monomial(width, 1, (0, 0, 0, 1))
    monom_X_4 = Monomial(width, -1, (0, 0, 0, 1))
    equation_H_4 = Equation(0, width, 1, [monom_H_4])
    equation_I_4 = Equation(1, width, 1, [monom_I_4])
    equation_T_4 = Equation(2, width, 1, [monom_T_4])
    equation_X_4 = Equation(3, width, 1, [monom_X_4])
    test_4 = Test(width, [equation_H_4, equation_I_4, equation_T_4, equation_X_4])
    return test_3, test_4, order

def hilln_benchmark_tests(repeat):
    """Perform a series of tests on hilln systems of ODEs (with both 3 and 4
    variables), for a given number of repetitions"""
    with open('hilln_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test_3, test_4, order = generate_hilln_tests()
            outfile.write(f'Order n: {order}\n')
            outfile.write('\nWith 3 variables:\n')
            outfile.write(str(test_3))
            test_3.run()
            outfile.write('\nSolution:\n')
            for laurent in test_3.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test_3.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test_3.all_substitutions)))
            outfile.write('\n\nWith 4 variables:\n')
            outfile.write(str(test_4))
            test_4.run()
            outfile.write('\nSolution:\n')
            for laurent in test_4.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test_4.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test_4.all_substitutions))) 
            outfile.write('\n\n')
    
def generate_selkov_test():
    """Create an instance of class Test that represents a selkov system 
    of ODEs, return it and the parameters a and b of the system"""
    width = 2 # two variables X(index 0) and Y(index 1)
    a = random.randint(-100, 100)
    b = random.randint(-100, 100)
    monomials_X = [Monomial(width, -1, (1, 0)), Monomial(width, a, (0, 1)), Monomial(width, 1, (2, 1))]
    monomials_Y = [Monomial(width, b, (0, 0)), Monomial(width, - a, (0, 1)), Monomial(width, -1, (2, 1))]
    eqX = Equation(0, width, 3, monomials_X)
    eqY = Equation(1, width, 3, monomials_Y)
    return Test(width, [eqX, eqY]), a, b

def selkov_benchmark_tests(repeat):
    """Perform a series of tests on selkov systems of ODEs, for a given number
    of repetitions"""
    with open('selkov_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test, a, b = generate_selkov_test()
            outfile.write(f'parameters: a = {a}, b = {b}\n')
            outfile.write(str(test))
            test.run()
            outfile.write('\nSoultion:\n')
            for laurent in test.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\n\n')
        
def generate_cubic_cycle_test():
    """Create an instance of class Test that represents a cubic cycle system 
    of ODEs"""
    width = random.randint(1, 15)
    equations = []
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 3
        monom = Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(Equation(i, width, 1, [monom]))
    return Test(width, equations)

def cubic_cycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic cycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_cycle_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_cubic_cycle_test()
            outfile.write(str(test))
            test.run()
            outfile.write('\nSoultion:\n')
            for laurent in test.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\n\n')
            
def generate_cubic_bicycle_test():
    """Create an instance of class Test that represents a cubic bicycle system 
    of ODEs"""
    width = random.randint(1, 10)
    equations = []
    l = [0] * width
    for i in range(width):
        l[i - 1] = 3
        monom_left = Monomial(width, 1, tuple(l))
        l[i - 1] = 0
        l[(i + 1) % width] = 3
        monom_right = Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(Equation(i, width, 2, [monom_left, monom_right]))
    return Test(width, equations)

def cubic_bicycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic bicycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_bicycle_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_cubic_bicycle_test()
            outfile.write(str(test))
            test.run()
            outfile.write('\nSoultion:\n')
            for laurent in test.optimal_solution:
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\n\n')
    
    

    
