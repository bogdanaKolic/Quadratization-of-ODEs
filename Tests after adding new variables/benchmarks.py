# -*- coding: utf-8 -*-

import quadratization
import random
import time
import itertools

def generate_random_circular_test():
    """Create an instance of class Test that represents a circular system 
    of ODEs"""
    width = random.randint(1, 15) # n
    degree = random.randint(1, 50) # k
    equations = []
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = degree
        monom = quadratization.Monomial(width, 1, tuple(l))
        equations.append(quadratization.Equation(i, width, 1, [monom]))
        l[(i + 1) % width] = 0
    return quadratization.QuadratizationProblem(width, equations)
    

def random_circular_benchmark_tests(repeat):
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
            
def generate_hardk_test(add_subs = False):
    """Create an instance of class Test that represents a hardk system 
    of ODEs, return it and the degree k of the system"""
    width = 3 # variables A(index 0), B(index 1) and C(index 2)
    degree = random.randint(1, 5) # k
    monomials_of_A = [quadratization.Monomial(width, 1, (0, 0, degree)), quadratization.Monomial(width, 1, (2, 2, degree - 1))]
    eqA = quadratization.Equation(0, width, 2, monomials_of_A)
    monomials_of_B = [quadratization.Monomial(width, 1, (2, 0, 0))]
    eqB = quadratization.Equation(1, width, 1, monomials_of_B)
    monomials_of_C = [quadratization.Monomial(width, 1, (0, 2, 0))]
    eqC = quadratization.Equation(2, width, 1, monomials_of_C)
    added_substitutions = set()
    if add_subs:
        added_substitutions = variables(width, [eqA, eqB, eqC], degree)
    return quadratization.QuadratizationProblem(width, [eqA, eqB, eqC], added_substitutions), degree

def random_hardk_benchmark_tests(repeat):
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
            
def generate_random_monomn_test():
    """Create an instance of class Test that represents a monomn system 
    of ODEs"""
    width = random.randint(1, 9)
    equations = []
    common_monomial = quadratization.Monomial(width, 1, (2, ) * width)
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 2
        monom = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 2, [common_monomial, monom]))
    return quadratization.QuadratizationProblem(width, equations)

def random_monomn_benchmark_tests(repeat):
    """Perform a series of tests on monomn systems of ODEs, for a given
    number of repetitions"""
    with open('monom_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_random_monomn_test()
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

def generate_random_hilln_tests():
    """Create two instances of class Test that represent hill systems
    of ODEs, one with 3, and one with 4 variables; return them and their order"""
    order = random.randint(1, 100) # n
    width = 3 # variables H(index 0), I(index 1) and T(index 2)
    monom_H_3 = quadratization.Monomial(width, order, (0, 2, order - 1))
    monom_I_3 = quadratization.Monomial(width, - order, (0, 2, order - 1))
    monom_T_3 = quadratization.Monomial(width, 1, (0, 0, 0))
    equation_H_3 = quadratization.Equation(0, width, 1, [monom_H_3])
    equation_I_3 = quadratization.Equation(1, width, 1, [monom_I_3])
    equation_T_3 = quadratization.Equation(2, width, 1, [monom_T_3])
    test_3 = quadratization.QuadratizationProblem(width, [equation_H_3, equation_I_3, equation_T_3])
    width = 4 # variables H(index 0), I(index 1) and T(index 2) and X(index 4)
    monom_H_4 = quadratization.Monomial(width, order, (0, 2, order - 1, 1))
    monom_I_4 = quadratization.Monomial(width, - order, (0, 2, order - 1, 1))
    monom_T_4 = quadratization.Monomial(width, 1, (0, 0, 0, 1))
    monom_X_4 = quadratization.Monomial(width, -1, (0, 0, 0, 1))
    equation_H_4 = quadratization.Equation(0, width, 1, [monom_H_4])
    equation_I_4 = quadratization.Equation(1, width, 1, [monom_I_4])
    equation_T_4 = quadratization.Equation(2, width, 1, [monom_T_4])
    equation_X_4 = quadratization.Equation(3, width, 1, [monom_X_4])
    test_4 = quadratization.QuadratizationProblem(width, [equation_H_4, equation_I_4, equation_T_4, equation_X_4])
    return test_3, test_4, order

def random_hilln_benchmark_tests(repeat):
    """Perform a series of tests on hilln systems of ODEs (with both 3 and 4
    variables), for a given number of repetitions"""
    with open('hilln_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test_3, test_4, order = generate_random_hilln_tests()
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
    
def generate_selkov_test(add_subs = False):
    """Create an instance of class Test that represents a selkov system 
    of ODEs, return it and the parameters a and b of the system"""
    width = 2 # two variables X(index 0) and Y(index 1)
    a = random.randint(-100, 100)
    b = random.randint(-100, 100)
    monomials_X = [quadratization.Monomial(width, -1, (1, 0)), quadratization.Monomial(width, a, (0, 1)), quadratization.Monomial(width, 1, (2, 1))]
    monomials_Y = [quadratization.Monomial(width, b, (0, 0)), quadratization.Monomial(width, - a, (0, 1)), quadratization.Monomial(width, -1, (2, 1))]
    eqX = quadratization.Equation(0, width, 3, monomials_X)
    eqY = quadratization.Equation(1, width, 3, monomials_Y)
    added_substitutions = set()
    if add_subs:
        added_substitutions = variables(width, [eqX, eqY], 2)
    return quadratization.QuadratizationProblem(width, [eqX, eqY], added_substitutions), a, b

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
        
def generate_random_cubic_cycle_test():
    """Create an instance of class Test that represents a cubic cycle system 
    of ODEs"""
    width = random.randint(1, 15)
    equations = []
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 3
        monom = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 1, [monom]))
    return quadratization.QuadratizationProblem(width, equations)

def random_cubic_cycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic cycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_cycle_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_random_cubic_cycle_test()
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
            
def generate_random_cubic_bicycle_test():
    """Create an instance of class Test that represents a cubic bicycle system 
    of ODEs"""
    width = random.randint(1, 20)
    equations = []
    l = [0] * width
    for i in range(width):
        l[i - 1] = 3
        monom_left = quadratization.Monomial(width, 1, tuple(l))
        l[i - 1] = 0
        l[(i + 1) % width] = 3
        monom_right = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 2, [monom_left, monom_right]))
    return quadratization.QuadratizationProblem(width, equations)

def random_cubic_bicycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic bicycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_bicycle_benchmark_tests.txt', 'a') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            test = generate_random_cubic_bicycle_test()
            outfile.write(str(test))
            test.run()
            outfile.write('\nSolution:\n')
            for i, laurent in enumerate(test.optimal_solution):
                outfile.write('y' + str(i) + ' = ')
                outfile.write(str(laurent.variables))
                outfile.write('\n')
            outfile.write('optimal number of substitutions: ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions: ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\n\n')
    
# following functions were used for testing the second algorihtm

def format_time(sec):
    mins = sec // 60
    sec %= 60
    h = mins // 60
    mins %= 60
    s = ''
    if h != 0:
        if h % 10 == 1 and h % 100 != 11:
            s += f'{h} hour '
        else:
            s += f'{h} hours '
    if mins != 0:
        if mins % 10 == 1 and mins % 100 != 11:
            s += f'{mins} minute '
        else:
            s += f'{mins} minutes '
    if sec % 10 == 1 and sec % 100 != 11:
            s += f'{sec} second '
    else:
        s += f'{sec} seconds '
    return s

def generate_circular_test(width, add_subs = False):
    """Create an instance of class Test that represents a circular system 
    of ODEs"""
    degree = random.randint(1, 10) # k
    equations = []
    l = [0]*width
    for i in range(width):
        l[(i + 1) % width] = degree
        monom = quadratization.Monomial(width, 1, tuple(l))
        equations.append(quadratization.Equation(i, width, 1, [monom]))
        l[(i + 1) % width] = 0
    added_substitutions = set()
    if add_subs:
        added_substitutions = variables(width, equations, degree)
    return quadratization.QuadratizationProblem(width, equations, added_substitutions), degree

def circular_benchmark_tests(repeat):
    """Perform a series of tests on circular systems of ODEs, for a given
    number of repetitions"""
    with open('circular_benchmark_tests_sub_from_derivatives_and_deg3.txt', 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1} :\n')
            t0 = time.time()
            test, degree = generate_circular_test(i, True)
            test.run()
            t1 = time.time()
            outfile.write('degree k = ')
            outfile.write(str(degree) + '\n')
            outfile.write('number of variables n = ' + str(test.width))
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test.optional_substitutions if y not in test.optimal_solution]):
                outfile.write(f'y{j + test.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')

def hardk_benchmark_tests(repeat):
    """Perform a series of tests on hardk systems of ODEs, for a given
    number of repetitions"""
    with open('hardk_benchmark_tests_sub_from_derivatives.txt', 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1} :\n')
            t0 = time.time()
            test, degree = generate_hardk_test(False)
            test.run()
            t1 = time.time()
            outfile.write('degree k = ')
            outfile.write(str(degree) + '\n')
            outfile.write('number of variables n = ' + str(test.width))
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test.optional_substitutions if y not in test.optimal_solution]):
                outfile.write(f'y{j + test.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')

def generate_monomn_test(width, add_subs = False):
    """Create an instance of class Test that represents a monomn system 
    of ODEs"""
    equations = []
    common_monomial = quadratization.Monomial(width, 1, (2, ) * width)
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 2
        monom = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 2, [common_monomial, monom]))
    added_substitutions = set()
    if add_subs:
        added_substitutions = variables(width, equations, 2)
    return quadratization.QuadratizationProblem(width, equations, added_substitutions)

def monomn_benchmark_tests(repeat):
    """Perform a series of tests on monomn systems of ODEs, for a given
    number of repetitions"""
    with open('monomn_benchmark_tests_sub_from_derivatives.txt', 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1} :\n')
            t0 = time.time()
            test = generate_monomn_test(i + 1, False)
            test.run()
            t1 = time.time()
            outfile.write('number of variables n = ' + str(test.width))
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test.optional_substitutions if y not in test.optimal_solution]):
                outfile.write(f'y{j + test.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')          

def generate_hilln_test(order, width, all_subs = False):
    """Create two instances of class Test that represent hill systems
    of ODEs, one with 3, and one with 4 variables; return them and their order"""
    # width = 3 variables H(index 0), I(index 1) and T(index 2)
    if width == 3:
        monom_H_3 = quadratization.Monomial(width, order, (0, 2, order - 1))
        monom_I_3 = quadratization.Monomial(width, - order, (0, 2, order - 1))
        monom_T_3 = quadratization.Monomial(width, 1, (0, 0, 0))
        equation_H_3 = quadratization.Equation(0, width, 1, [monom_H_3])
        equation_I_3 = quadratization.Equation(1, width, 1, [monom_I_3])
        equation_T_3 = quadratization.Equation(2, width, 1, [monom_T_3])
        added_substitutions = set()
        if all_subs:
            added_substitutions = variables(width, [equation_H_3, equation_I_3, equation_T_3], order)
        return quadratization.QuadratizationProblem(width, [equation_H_3, equation_I_3, equation_T_3], added_substitutions)
    # width = 4  variables H(index 0), I(index 1) and T(index 2) and X(index 4)
    monom_H_4 = quadratization.Monomial(width, order, (0, 2, order - 1, 1))
    monom_I_4 = quadratization.Monomial(width, - order, (0, 2, order - 1, 1))
    monom_T_4 = quadratization.Monomial(width, 1, (0, 0, 0, 1))
    monom_X_4 = quadratization.Monomial(width, -1, (0, 0, 0, 1))
    equation_H_4 = quadratization.Equation(0, width, 1, [monom_H_4])
    equation_I_4 = quadratization.Equation(1, width, 1, [monom_I_4])
    equation_T_4 = quadratization.Equation(2, width, 1, [monom_T_4])
    equation_X_4 = quadratization.Equation(3, width, 1, [monom_X_4])
    added_substitutions = set()
    if all_subs:
        added_substitutions = variables(width, [equation_H_4, equation_I_4, equation_T_4, equation_X_4], order)
    return quadratization.QuadratizationProblem(width, [equation_H_4, equation_I_4, equation_T_4, equation_X_4], added_substitutions)

def hilln_benchmark_tests(repeat):
    """Perform a series of tests on hilln systems of ODEs (with both 3 and 4
    variables), for a given number of repetitions"""
    with open('hilln_benchmark_tests_sub_from_derivatives.txt', 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            #order = random.randint(1, 10)
            order = i + 1
            outfile.write(f'Order n: {order}\n')
            outfile.write('\nWith 3 variables:\n')
            t0 = time.time()
            test_3 = generate_hilln_test(order, 3)
            test_3.run()
            t1 = time.time()
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test_3.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test_3.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test_3))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test_3.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test_3.optional_substitutions if y not in test_3.optimal_solution]):
                outfile.write(f'y{j + test_3.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')
            outfile.write('\n\nWith 4 variables:\n')
            t0 = time.time()
            test_4 = generate_hilln_test(order, 4)
            test_4.run()
            t1 = time.time()
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test_4.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test_4.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test_4))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test_4.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test_4.optional_substitutions if y not in test_4.optimal_solution]):
                outfile.write(f'y{j + test_4.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')
            
def selkov_benchmark_tests2(repeat):
    """Perform a series of tests on selkov systems of ODEs, for a given number
    of repetitions"""
    with open('selkov_benchmark_tests_sub_from_derivatives.txt', 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            t0 = time.time()
            test, a, b = generate_selkov_test()
            test.run()
            t1 = time.time()
            outfile.write(f'parameters: a = {a}, b = {b}\n')
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test.optional_substitutions if y not in test.optimal_solution]):
                outfile.write(f'y{j + test.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')  
            
def generate_cubic_cycle_test(width, add_subs = False):
    """Create an instance of class Test that represents a cubic cycle system 
    of ODEs"""
    equations = []
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 3
        monom = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 1, [monom]))
    added_substitutions = set()
    if add_subs:
        added_substitutions = variables(width, equations, 3)
    return quadratization.QuadratizationProblem(width, equations, added_substitutions)

def cubic_cycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic cycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_cycle_benchmark_tests_sub_from_derivatives.txt', 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            t0 = time.time()
            test = generate_cubic_cycle_test(i)
            test.run()
            t1 = time.time()
            outfile.write('number of variables n = ' + str(test.width))
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test.optional_substitutions if y not in test.optimal_solution]):
                outfile.write(f'y{j + test.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')
            
def generate_cubic_bicycle_test(width, add_subs = False):
    """Create an instance of class Test that represents a cubic bicycle system 
    of ODEs"""
    equations = []
    l = [0] * width
    for i in range(width):
        l[i - 1] = 3
        monom_left = quadratization.Monomial(width, 1, tuple(l))
        l[i - 1] = 0
        l[(i + 1) % width] = 3
        monom_right = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 2, [monom_left, monom_right]))
    added_substitutions = set()
    if add_subs:
        added_substitutions = variables(width, equations, 3)
    return quadratization.QuadratizationProblem(width, equations, added_substitutions)

def cubic_bicycle_benchmark_tests(repeat):
    """Perform a series of tests on cubic bicycle systems of ODEs, for a given 
    number of repetitions"""
    with open('cubic_bicycle_benchmark_tests_sub_from_derivatives.txt', 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            t0 = time.time()
            test = generate_cubic_bicycle_test(i)
            test.run()
            t1 = time.time()
            outfile.write('number of variables n = ' + str(test.width))
            outfile.write('\noptimal number of substitutions : ')
            outfile.write(str(test.min_length))
            outfile.write('\nnumber of all substitutions : ')
            outfile.write(str(len(test.all_substitutions))) 
            outfile.write('\ntime of execution : ')
            outfile.write(format_time(t1 - t0))
            outfile.write('\n' + str(test))
            outfile.write('\nSolution :\n')
            for j, laurent in enumerate(test.optimal_solution):
                outfile.write(f'y{j} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\nDisregarded substitutions :\n')
            for j, laurent in enumerate([y for y in test.optional_substitutions if y not in test.optimal_solution]):
                outfile.write(f'y{j + test.min_length} = ')
                outfile.write(str(laurent))
                outfile.write('\n')
            outfile.write('\n')
            
def variables(n, equations, degree):
    """ Function for creating additional substitution variables"""
    substitutions = set()
    l = [j for j in itertools.combinations_with_replacement([i for i in range(-3, 4)], n)]
    for comb in l: 
        for s in itertools.permutations(comb):
            substitutions.add(quadratization.Substitution(s, equations))
    return substitutions

