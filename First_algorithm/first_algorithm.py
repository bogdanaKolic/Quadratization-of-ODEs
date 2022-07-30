# -*- coding: utf-8 -*-

import random

max_coeff = 100 # maximal generated coefficient of a monomial
max_power = 30 # maximal generated power of a variable in a monomial
max_width = 3 # maximal number of original variables
max_length = 7 # maximal number of monomials a generated equation can have

class Monomial:
    """ Monomials are represented with a coefficient, number of variables, 
    and a tuple of the powers of each variable """
    def __init__(self, n = 0, c = 1, var = tuple()):
        self.width = n # number of variables
        self.coefficient = c
        self.variables = var
            
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
    
    def __hash__(self):
        return hash(self.variables)
    
    def __str__(self):
        return str(self.coefficient) + '*' + str(self.variables)      
    
class Laurent(Monomial):
    """Laurent monomials in addition have the index of the variable by which 
    they are divided and a pre-calculated derivative"""
    def __init__(self, monomial, index, equations):
        Monomial.__init__(self, monomial.width, monomial.coefficient)
        for i, deg in enumerate(monomial.variables):
            if i == index:
                self.variables += (monomial.variables[i] - 1, )
            else:
                self.variables += (monomial.variables[i], )
        self.derivative = self.calculate_derivative(equations)
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
        # s += '\nderivative: ' + str(self.derivative) + '\nindex: ' + str(self.index)
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
        self.monomials = monomials
    
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
        for e in self.equations:
            self.all_substitutions += e.calculate_substitutions(self.equations)
        self.current_substitutions = [] # list of substitutions currently cosidered as a quadratization
        self.optimal_solution = [y for y in self.all_substitutions]
        self.min_length = len(self.optimal_solution)
    
    
    def random_test(self):
        self.width = random.randint(1, max_width)
        self.equations = []
        for i in range(self.width):
            eq = Equation(i, self.width)
            eq.set_random(self.width)
            self.equations.append(eq)
        self.all_substitutions = [] 
        for e in self.equations:
            self.all_substitutions += e.calculate_substitutions(self.equations)
        self.current_substitutions = [] # list of substitutions currently cosidered as a quadratization
        self.optimal_solution = [y for y in self.all_substitutions]
        self.min_length = len(self.optimal_solution)
    
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
        self.all_substitutions = [] 
        for e in self.equations:
            self.all_substitutions += e.calculate_substitutions(self.equations)
        self.current_substitutions = [] # list of substitutions currently cosidered as a quadratization
        self.optimal_solution = [y for y in self.all_substitutions]
        self.min_length = len(self.optimal_solution)
        
    def is_quadratization(self):
        """ Checks if the set of current_substitutions works as a quadratization """
        search_space = set()
        # add 'a constant'
        search_space.add((0, ) * self.width)
        # add variables xi, xixj and xi^2
        l = [0] * self.width
        for i in range(self.width):
            l[i] += 1
            search_space.add(tuple(l))
            for j in range(self.width):
                l[j] += 1
                search_space.add(tuple(l))
                l[j] -= 1
            l[i] -= 1
        # add substitution variables and their products
        for first in self.current_substitutions:
            search_space.add(first.variables)
            for second in self.current_substitutions:
                search_space.add((Monomial.multiply(first, second)).variables)
        # add original variables multiplied by Laurent monomials
        for monom in self.current_substitutions:
            for i in range(self.width):
                search_space.add((monom.multiply_x(i)).variables)
        # compute the list of all monomials in the derivateives
        derivatives = []
        for poly in self.equations:
            for mono in poly.monomials:
                derivatives.append(mono.variables)
        for substitution in self.current_substitutions:
            for monom in substitution.derivative.monomials:
                derivatives.append(monom.variables)
        for derivative in derivatives:
            if derivative not in search_space:
                return False
        return True
    

    def reduce(self, position):
        """ Find out how many Laurent monomials can be neglected by using 
        a recursion on a set of all substitutions and considering the monomial 
        at a given position"""
        if position == len(self.all_substitutions):
            if len(self.current_substitutions) < self.min_length and self.is_quadratization():
                self.min_length = len(self.current_substitutions)
                self.optimal_solution = [y for y in self.current_substitutions]
        elif len(self.current_substitutions) >= self.min_length:
            return None
        else:
            self.reduce(position + 1)
            self.current_substitutions.append(self.all_substitutions[position])
            self.reduce(position + 1)
            self.current_substitutions.remove(self.all_substitutions[position])
    
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
        # s += 'original differential equations:\n'
        # for e in self.equations:
        #     s += str(e) + '\n'
        s += 'substitutions:\n'
        for a in self.all_substitutions:
            s += str(a) + '\n'
        return s
        
def main_random():
    """ Simulate quadratizing a set of equations """
    for _ in range(10):
        t = Test()
        t.random_test()
        t.run()
    

