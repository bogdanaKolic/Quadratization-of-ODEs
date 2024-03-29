# -*- coding: utf-8 -*-

import quadratization
import random
import time
import os

strategies = ['no additional substitutions', 'constant degree', 'roots of monom derivatives',
              'random from derivatives', 'from derivatives']

NUM_RUNS_BENCHMARKING = 10
NUM_RUNS_TIME_MEASUREMENT = 5
MAX_WIDTH_BENCHMARKS = 15
MAX_DEGREE_BENCHMARKS = 100
SELKOV_PARAMETERS_RANGE = 100

def format_time(sec):
    """Transform the number of seconds given into a string representation
    of time stating the number of hours, minutes and seconds"""
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

def write_to_file(iteration, test, time, outfile, degree = None, additional_string = None):
    """ A function for writing the results of a single test into a file"""
    outfile.write(f'Test {iteration} :\n')
    if not additional_string is None:
        outfile.write(additional_string)
    outfile.write('number of variables n = ' + str(test.width))
    if not degree is None:
        outfile.write('\ndegree k = ' + str(degree))
    outfile.write('\noptimal number of substitutions : ')
    outfile.write(str(test.min_length))
    outfile.write('\nnumber of all substitutions : ')
    outfile.write(str(len(test.all_substitutions))) 
    outfile.write('\ntime of execution : ')
    outfile.write(format_time(time))
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

def generate_circular_test(width, degree):
    """Create an instance of class QuadratizationProblem that represents a 
    circular system of ODEs, with given parameters or randomly; return the 
    object and the degree of the system"""
    if width is None:
       width = random.randint(1, MAX_WIDTH_BENCHMARKS)
    if degree is None:
        degree = random.randint(1, MAX_DEGREE_BENCHMARKS)
    equations = []
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = degree
        monom = quadratization.Monomial(width, 1, tuple(l))
        equations.append(quadratization.Equation(i, width, 1, [monom]))
        l[(i + 1) % width] = 0
    return quadratization.QuadratizationProblem(width, equations), degree       
            
def generate_hardk_test(width, degree):
    """Create an instance of class QuadratizationProblem that represents a 
    hardk system of ODEs, either of the given degree or randomly; return it 
    and the degree k of the system"""
    if degree is None:
       degree = random.randint(1, MAX_DEGREE_BENCHMARKS)
    width = 3 # variables A(index 0), B(index 1) and C(index 2)
    # degree = random.randint(1, 5) # k
    monomials_of_A = [quadratization.Monomial(width, 1, (0, 0, degree)), quadratization.Monomial(width, 1, (2, 2, degree - 1))]
    eqA = quadratization.Equation(0, width, 2, monomials_of_A)
    monomials_of_B = [quadratization.Monomial(width, 1, (2, 0, 0))]
    eqB = quadratization.Equation(1, width, 1, monomials_of_B)
    monomials_of_C = [quadratization.Monomial(width, 1, (0, 2, 0))]
    eqC = quadratization.Equation(2, width, 1, monomials_of_C)
    return quadratization.QuadratizationProblem(width, [eqA, eqB, eqC]), degree
            
def generate_monomn_test(width, degree = None):
    """Create an instance of class QuadratizationProblem that represents a 
    monomn system of ODEs, with the given or random number of equations, and
    return it"""
    if width is None:
        width = random.randint(1, MAX_WIDTH_BENCHMARKS)
    equations = []
    common_monomial = quadratization.Monomial(width, 1, (2, ) * width)
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 2
        monom = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 2, [common_monomial, monom]))
    return quadratization.QuadratizationProblem(width, equations), None

def generate_hilln_test(width, order):
    """Create an instance of class QuadratizationProblem that represents hill 
    systems of ODEs with the given number of equations, and of either given or
    randomly chosen order, then return both the object and the order"""
    if order is None:
        order = random.randint(1, MAX_DEGREE_BENCHMARKS)
    if width == 3: # variables H(index 0), I(index 1) and T(index 2)
        monom_H_3 = quadratization.Monomial(width, order, (0, 2, order - 1))
        monom_I_3 = quadratization.Monomial(width, - order, (0, 2, order - 1))
        monom_T_3 = quadratization.Monomial(width, 1, (0, 0, 0))
        equation_H_3 = quadratization.Equation(0, width, 1, [monom_H_3])
        equation_I_3 = quadratization.Equation(1, width, 1, [monom_I_3])
        equation_T_3 = quadratization.Equation(2, width, 1, [monom_T_3])
        return quadratization.QuadratizationProblem(width, [equation_H_3, equation_I_3, equation_T_3]), order
    # width = 4: variables H(index 0), I(index 1) and T(index 2) and X(index 4)
    monom_H_4 = quadratization.Monomial(width, order, (0, 2, order - 1, 1))
    monom_I_4 = quadratization.Monomial(width, - order, (0, 2, order - 1, 1))
    monom_T_4 = quadratization.Monomial(width, 1, (0, 0, 0, 1))
    monom_X_4 = quadratization.Monomial(width, -1, (0, 0, 0, 1))
    equation_H_4 = quadratization.Equation(0, width, 1, [monom_H_4])
    equation_I_4 = quadratization.Equation(1, width, 1, [monom_I_4])
    equation_T_4 = quadratization.Equation(2, width, 1, [monom_T_4])
    equation_X_4 = quadratization.Equation(3, width, 1, [monom_X_4])
    return quadratization.QuadratizationProblem(width, [equation_H_4, equation_I_4, equation_T_4, equation_X_4]), order
    
def generate_selkov_test():
    """Create an instance of class QuadratizationProblem that represents a 
    selkov system of ODEs, return it and the randomly chosen parameters
    a and b of the system"""
    width = 2 # two variables X(index 0) and Y(index 1)
    a = random.randint(-SELKOV_PARAMETERS_RANGE, SELKOV_PARAMETERS_RANGE)
    b = random.randint(-SELKOV_PARAMETERS_RANGE, SELKOV_PARAMETERS_RANGE)
    monomials_X = [quadratization.Monomial(width, -1, (1, 0)), quadratization.Monomial(width, a, (0, 1)), quadratization.Monomial(width, 1, (2, 1))]
    monomials_Y = [quadratization.Monomial(width, b, (0, 0)), quadratization.Monomial(width, - a, (0, 1)), quadratization.Monomial(width, -1, (2, 1))]
    eqX = quadratization.Equation(0, width, 3, monomials_X)
    eqY = quadratization.Equation(1, width, 3, monomials_Y)
    return quadratization.QuadratizationProblem(width, [eqX, eqY]), a, b
        
def generate_cubic_cycle_test(width, degree = None):
    """Create and return an instance of class QuadratizationProblem that 
    represents a cubic cycle system of ODEs, with either given or randomly chosen
    number of equations"""
    if width is None:
        width = random.randint(1, MAX_WIDTH_BENCHMARKS)
    equations = []
    l = [0] * width
    for i in range(width):
        l[(i + 1) % width] = 3
        monom = quadratization.Monomial(width, 1, tuple(l))
        l[(i + 1) % width] = 0
        equations.append(quadratization.Equation(i, width, 1, [monom]))
    return quadratization.QuadratizationProblem(width, equations), None

def generate_cubic_bicycle_test(width, degree = None):
    """Create and return an instance of class QuadratizationProblem that 
    represents a cubic bicycle system of ODEs, with either given or randomly 
    chosen number of equations"""
    if width is None:
        width = random.randint(1, MAX_WIDTH_BENCHMARKS)
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
    return quadratization.QuadratizationProblem(width, equations), None

    
# following functions were used for testing the second algorihtm

def define_outfile_name(test, strategy, width = None, degree = None, degree_sub = 0):
    """ Based on the passed arguments, determine the name of the output file
    for a test"""
    str_to_name = {'from derivatives' : '_sub_from_derivatives',
                     'random from derivatives' : '_random_sub_from_derivatives',
                     'constant degree' : f'_degree_{degree_sub}_substitutions',
                     'roots of monom derivatives' : '_derivative_roots',
                     'optimal quadratization' : '_with_optquad'}
    if width != None and degree != None:
        test += f'({width},{degree})'
    elif width != None:
        test += f'({width})'
    elif degree != None:
        test += f'_degree{degree}'
    if strategy in str_to_name:
        test += str_to_name[strategy]
    test += '.txt'
    return test

def benchmark_tests(benchmark, repeat, strategy = 'no additional substitutions,', width = None, degree = None, degree_sub = 0, optquad = None):
    """Perform a series of tests on systems of ODEs, with given parameters"""
    benchmark_to_function = {'circular' : generate_circular_test,
                             'hard' : generate_hardk_test,
                             'monom' : generate_monomn_test,
                             'cubic_cycle' : generate_cubic_cycle_test,
                             'cubic_bicycle' : generate_cubic_bicycle_test}
    filename = define_outfile_name(benchmark, strategy, width, degree, degree_sub = degree_sub)
    with open(filename, 'w') as outfile:
        for i in range(repeat):
            t0 = time.time()
            test, degree = benchmark_to_function[benchmark](width, degree)
            test.run(strategy, degree_sub, optquad)
            t1 = time.time()
            write_to_file(i + 1, test, t1 - t0, outfile, degree)       
            
def hilln_benchmark_tests(repeat, strategy = 'no additional substitutions,', order = None, degree_sub = 0, optquad = None):
    """Perform a series of tests on hilln systems of ODEs (with both 3 and 4
    variables), with given parameters"""
    filename = define_outfile_name('hill', strategy, degree = order, degree_sub = degree_sub)
    with open(filename, 'w') as outfile:
        for i in range(repeat):
            outfile.write(f'Test {i + 1}:\n')
            outfile.write()
            t0 = time.time()
            test_3 = generate_hilln_test(3, order)
            test_3.run(strategy, degree_sub, optquad)
            t1 = time.time()
            write_to_file(i + 1, test_3, t1 - t0, outfile, additional_string = 'With 3 variables:\n')
            t0 = time.time()
            test_4 = generate_hilln_test(4, order)
            test_4.run(strategy, degree_sub, optquad)
            t1 = time.time()
            write_to_file(i + 1, test_4, t1 - t0, outfile, additional_string = 'With 4 variables:\n')
                
def selkov_benchmark_tests(repeat, strategy = 'no additional substitutions,', degree_sub = 0, optquad = None):
    """Perform a series of tests on selkov systems of ODEs, with given parameters"""
    filename = define_outfile_name('selkov', strategy, degree_sub = degree_sub)
    with open(filename, 'w') as outfile:
        for i in range(repeat):
            t0 = time.time()
            test, a, b = generate_selkov_test()
            test.run(strategy, degree_sub, optquad)
            t1 = time.time()
            write_to_file(i + 1, test, t1 - t0, outfile, additional_string = f'parameters: a = {a}, b = {b}\n') 

def hiv_benchmark_test(repeat, strategy = 'no additional substitutions,', degree_sub = 0, optquad = None):
    """Perform a  series of tests on a hiv system of ODEs, with given parameters"""
    filename = define_outfile_name('hiv', strategy, degree_sub = degree_sub)
    with open(filename, 'w') as outfile:
        for i in range(repeat):
            t0 = time.time()
            test = quadratization.QuadratizationProblem()
            test.load_from_file('input_files/hiv.txt')
            test.run(strategy, degree_sub, optquad)
            t1 = time.time()
            write_to_file(i + 1, test, t1 - t0, outfile)

def compare_must_have_sub(outfile, t_without, t_with, degree_sub):
    """Perform a single comparison regarding the number of must-have substitution
    variables"""
    t_without.precompute_part2()
    t_without.remove_impossible_substitutions()
    t_without.find_must_have_substitutions()
    t_with.variables_constant_degree(degree_sub)
    t_with.precompute_part2()
    t_with.remove_impossible_substitutions()
    t_with.find_must_have_substitutions()
    outfile.write(f'number of variables: {t_with.width}\n\n')
    outfile.write(f'number of must-have substitutions without addition: {len(t_without.all_substitutions) - len(t_without.optional_substitutions)}\n')
    outfile.write(f'number of must-have substitutions with addition: {len(t_with.all_substitutions) - len(t_with.optional_substitutions)}\n\n')
        

def count_must_have_sub(degree_sub):
    """Perform a series of tests to compare the number of must-have 
    substitutions with and without additional variables in the system"""
    with open(f'count_must-have-substitutions_degree_{degree_sub}.txt', 'w') as outfile:
        outfile.write('Circular\n')
        for i in range(NUM_RUNS_BENCHMARKING):
            outfile.write(f'Test {i} :\n')
            degree = random.randint(1, MAX_DEGREE_BENCHMARKS)
            outfile.write(f'degree k: {degree}\n')
            t_without, _ = generate_circular_test(i + 1, degree)
            t_with, _ = generate_circular_test(i + 1, degree)
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)    
        outfile.write('Hardk\n')
        for i in range(NUM_RUNS_BENCHMARKING):
            outfile.write(f'Test {i + 1} :\n')
            outfile.write(f'degree k: {i + 1}\n')
            t_without, _ = generate_hardk_test(i + 1)
            t_with, _ = generate_hardk_test(i + 1)
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)
        outfile.write('Monomn\n')
        for i in range(NUM_RUNS_BENCHMARKING):
            outfile.write(f'Test {i + 1} :\n')
            t_without, _ = generate_monomn_test(i + 1)
            t_with, _ = generate_monomn_test(i + 1)
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)
        outfile.write('Selkov\n')
        for i in range(NUM_RUNS_BENCHMARKING):
            outfile.write(f'Test {i + 1}:\n')
            t_without, _, _ = generate_selkov_test()
            t_with, _, _ = generate_selkov_test()
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)
        outfile.write('Cubic cycle\n')
        for i in range(NUM_RUNS_BENCHMARKING):
            outfile.write(f'Test {i + 1} :\n')
            t_without, _ = generate_cubic_cycle_test(i)
            t_with, _ = generate_cubic_cycle_test(i)
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)
        outfile.write('Cubic bicycle\n')
        for i in range(NUM_RUNS_BENCHMARKING):
            outfile.write(f'Test {i + 1} :\n')
            t_without, _ = generate_cubic_bicycle_test(i)
            t_with, _ = generate_cubic_bicycle_test(i)
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)
        outfile.write('Hiv\n')
        t_without = quadratization.QuadratizationProblem()
        t_without.load_from_file('input_files/hiv.txt')
        t_with = quadratization.QuadratizationProblem()
        t_with.load_from_file('input_files/hiv.txt')
        compare_must_have_sub(outfile, t_without, t_with, degree_sub)
        outfile.write('Hilln\n')
        for i in range(NUM_RUNS_BENCHMARKING):
            outfile.write(f'Test {i + 1} :\n')
            outfile.write(f'order: {i}\n')
            t_without, _ = generate_hilln_test(3, i)
            t_with, _ = generate_hilln_test(3, i)
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)
            t_without, _= generate_hilln_test(4, i)
            t_with, _= generate_hilln_test(4, i)
            compare_must_have_sub(outfile, t_without, t_with, degree_sub)
            
def test_from_file(foldername, system, strategy = 'no additional variables', degree_sub = 0, optquad = None, iteration = 1):
    """Perform a test on the data loaded from a file and return the execution time"""
    outfilename = define_outfile_name(system[: -4], strategy, degree_sub = degree_sub)
    with open(outfilename, 'a') as outfile:
        t0 = time.time()
        test = quadratization.QuadratizationProblem()
        test.load_from_file(foldername + '/' + system)
        test.run(strategy, degree_sub, optquad)
        t1 = time.time()
        write_to_file(iteration, test, t1 - t0, outfile)
    return t1 - t0

def tests_with_average_time(foldername, strategy = 'no additional substitutions'):
    """Perform tests on the set of benchmark problems stored in files, measure the 
    average time and the absolute error, given the strategy and the path to the files"""
    for file in os.listdir(foldername):
        t = 0
        times = []
        for i in range(NUM_RUNS_TIME_MEASUREMENT):
            dt = test_from_file(foldername, file, strategy, iteration = i + 1)
            t += dt
            times.append(dt)
        outfilename = define_outfile_name(file[: - 4], strategy)
        with open(outfilename, 'a') as outfile:
            outfile.write(f'average time: {format_time(t/NUM_RUNS_TIME_MEASUREMENT)}\n' )
            differences = [abs(x - t/NUM_RUNS_TIME_MEASUREMENT) for x in times]
            outfile.write(f'absolute error: {format_time(max(differences))}\n')

