# -*- coding: utf-8 -*-

import random
import first_algorithm as falg

def main_from_file(filename):
    """ Perform a test on the system stored in a file"""
    t = falg.Test()
    t.load_from_file(filename)
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
        monom = falg.Monomial(width, 1, tuple(l))
        equations.append(falg.Equation(i, width, 1, [monom]))
        l[(i + 1) % width] = 0
    return falg.Test(width, equations)
    

def circular_benchmark_tests(repeat):
    """Perform a series of tests on circular systems of ODEs, for a given
    number of repetitions"""
    with open('circular_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_circular_test()
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
            
def generate_hardk_test():
    """Create an instance of class Test that represents a hardk system 
    of ODEs, return it and the degree k of the system"""
    width = 3 # variables A(index 0), B(index 1) and C(index 2)
    degree = random.randint(1, 50) # k
    monomials_of_A = [falg.Monomial(width, 1, (0, 0, degree)), falg.Monomial(width, 1, (2, 2, degree - 1))]
    eqA = falg.Equation(0, width, 2, monomials_of_A)
    monomials_of_B = [falg.Monomial(width, 1, (2, 0, 0))]
    eqB = falg.Equation(1, width, 1, monomials_of_B)
    monomials_of_C = [falg.Monomial(width, 1, (0, 2, 0))]
    eqC = falg.Equation(2, width, 1, monomials_of_C)
    return falg.Test(width, [eqA, eqB, eqC]), degree

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
    common_monomial = falg.Monomial(width, 1, (2, ) * width)
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 2
        monom = falg.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(falg.Equation(i, width, 2, [common_monomial, monom]))
    return falg.Test(width, equations)

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
    monom_H_3 = falg.Monomial(width, order, (0, 2, order - 1))
    monom_I_3 = falg.Monomial(width, - order, (0, 2, order - 1))
    monom_T_3 = falg.Monomial(width, 1, (0, 0, 0))
    equation_H_3 = falg.Equation(0, width, 1, [monom_H_3])
    equation_I_3 = falg.Equation(1, width, 1, [monom_I_3])
    equation_T_3 = falg.Equation(2, width, 1, [monom_T_3])
    test_3 = falg.Test(width, [equation_H_3, equation_I_3, equation_T_3])
    width = 4 # variables H(index 0), I(index 1) and T(index 2) and X(index 4)
    monom_H_4 = falg.Monomial(width, order, (0, 2, order - 1, 1))
    monom_I_4 = falg.Monomial(width, - order, (0, 2, order - 1, 1))
    monom_T_4 = falg.Monomial(width, 1, (0, 0, 0, 1))
    monom_X_4 = falg.Monomial(width, -1, (0, 0, 0, 1))
    equation_H_4 = falg.Equation(0, width, 1, [monom_H_4])
    equation_I_4 = falg.Equation(1, width, 1, [monom_I_4])
    equation_T_4 = falg.Equation(2, width, 1, [monom_T_4])
    equation_X_4 = falg.Equation(3, width, 1, [monom_X_4])
    test_4 = falg.Test(width, [equation_H_4, equation_I_4, equation_T_4, equation_X_4])
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
    monomials_X = [falg.Monomial(width, -1, (1, 0)), falg.Monomial(width, a, (0, 1)), falg.Monomial(width, 1, (2, 1))]
    monomials_Y = [falg.Monomial(width, b, (0, 0)), falg.Monomial(width, - a, (0, 1)), falg.Monomial(width, -1, (2, 1))]
    eqX = falg.Equation(0, width, 3, monomials_X)
    eqY = falg.Equation(1, width, 3, monomials_Y)
    return falg.Test(width, [eqX, eqY]), a, b

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
            outfile.write('\nSolution:\n')
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
        monom = falg.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(falg.Equation(i, width, 1, [monom]))
    return falg.Test(width, equations)

def cubic_cycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic cycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_cycle_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_cubic_cycle_test()
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
            
def generate_cubic_bicycle_test():
    """Create an instance of class Test that represents a cubic bicycle system 
    of ODEs"""
    width = random.randint(1, 10)
    equations = []
    l = [0] * width
    for i in range(width):
        l[i - 1] = 3
        monom_left = falg.Monomial(width, 1, tuple(l))
        l[i - 1] = 0
        l[(i + 1) % width] = 3
        monom_right = falg.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(falg.Equation(i, width, 2, [monom_left, monom_right]))
    return falg.Test(width, equations)

def cubic_bicycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic bicycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_bicycle_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_cubic_bicycle_test()
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

