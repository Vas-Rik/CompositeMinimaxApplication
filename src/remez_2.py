from itertools import combinations
from matplotlib import pyplot as plt
from numpy.polynomial import Chebyshev
import sympy as sp
from scipy.optimize import root_scalar
from mpmath import mpf, mpc, mp
import numpy as np
from decimal import Decimal, getcontext
import mpmath as mp
import sys
import re


# Set precision to handle very small ranges
mp.mp.dps = 7000
getcontext().prec = 7000


def generate_chebyshev_points(n, lower, upper):
    """
    Generates a set of Chebyshev points spaced in the range [lower, upper].
    :param n: number of points
    :param lower: lower limit (Decimal)
    :param upper: upper limit (Decimal)
    :return: a list of Chebyshev points that are in the range [lower, upper]
    """
    index = np.arange(1, n + 1)
    range_ = abs(upper - lower)
    return [(mpf(0.5) * (mpf(mp.cos((2 * i - 1) / (2 * n) * mp.pi)) + mpf(1))) * range_ + lower for i in
            index]


def generate_points_in_ranges(ranges, total_num_points):
    """
    Generate Chebyshev points within given ranges to match the total number of points.
    :param ranges: A list of tuples that hold the ranges over which points have to be generated
    :param total_num_points: The total number of points that should be generated
    :return: A numpy array of Chebyshev points
    """
    total_length = sum(mpf(b) - mpf(a) for a, b in ranges)
    points = []

    for a, b in ranges:
        length = mpf(b) - mpf(a)
        num_points = round(total_num_points * (length / total_length))

        if num_points > 0:
            range_points = generate_chebyshev_points(num_points, mpf(a), mpf(b))
            points.extend(range_points)

    # Adjust the number of points to match total_num_points exactly
    if len(points) != total_num_points:
        if len(points) > total_num_points:
            points = points[:total_num_points]
        else:
            while len(points) < total_num_points:
                for a, b in ranges:
                    if len(points) >= total_num_points:
                        break
                    extra_points = generate_chebyshev_points(1, mpf(a), mpf(b))
                    if extra_points:
                        points.append(extra_points[0])

    points.sort()
    return np.array([p for p in points], dtype=mpf)


# Select n+1 points with alternating signs and maximum sum
def find_optimal_points(sets_of_extremes, r_i_diff):
    """
    Find the combination n+1 of points with alternating signs and maximum sum.
    :param sets_of_extremes: the set of total extreme points
    :param r_i_diff: the function for which the points have to be chosen
    """
    max_sum = -np.inf
    best_combination = None

    for extremes in sets_of_extremes:
        diffs = [r_i_diff(extreme) for extreme in extremes]
        if all(d1 * d2 < 0 for d1, d2 in zip(diffs, diffs[1:])):
            current_sum = sum(diffs)
            if current_sum > max_sum:
                max_sum = current_sum
                best_combination = extremes

    return best_combination


def find_combinations(X, n):
    """
    Generate all combinations of n+1 points from an ordered set X.

    :param X: The ordered list of points.
    :param n: The number of points in each combination minus one.

    Returns:
    list of tuples: Each tuple is a combination of n+1 points.
    """
    # Generate combinations of indices
    comb_indices = combinations(range(len(X)), n + 1)

    # Create combinations of points based on the indices
    comb_points = [tuple(X[i] for i in comb) for comb in comb_indices]

    return comb_points


def find_zeros(first_derivative, a, b):
    """
    Find zeros of the derivative of a Chebyshev polynomial in the interval [a, b].
    """
    x_sym = sp.symbols('x', real=True)

    # p_deriv_func = sp.lambdify(x_sym, first_derivative, 'numpy')
    p_deriv_func = first_derivative

    results = []
    zeros = []
    # Search for roots in the interval [a, b] by sampling
    num_samples = 3200
    # show_graph = False
    sample_points = mp.linspace(a, b, num_samples)
    # x = np.linspace(a, b, 50)
    #
    # if show_graph:
    #     y_exact = np.array([p_deriv_func(x_i) for x_i in x])
    #     plt.plot(x, y_exact, label='derivative function')
    #     plt.legend()
    #     plt.xlim(float(a), float(-0.9999999999999998))
    #     plt.title(f'Derivative temp')
    #     plt.show()

    # print("check_point_3.1")
    for i in range(num_samples - 1):
        x1, x2 = sample_points[i], sample_points[i + 1]
        if p_deriv_func(x1) * p_deriv_func(x2) < 0:
            # Use scipy's root_scalar to find the root within this interval
            try:
                # result = root_scalar(p_deriv_func, bracket=[x1, x2])
                zero = mp.findroot(p_deriv_func, (x1, x2), solver='anderson')
                zeros.append(zero)
            except ValueError:
                print("oops")

            # if result.converged:
            #     results.append(mpf(result.root))

    return zeros


def sign(x):
    """Returns the sign of x."""
    return (x >= 0) - (x < 0)


# ChooseNewNodes implements Algorithm 3 of High-Precision Bootstrapping
# of RNS-CKKS Homomorphic Encryption Using Optimal Minimax Polynomial
# But instead of concavity I use the sign function here
def ChooseNewXpoints(B, d, r):
    """
    Returns d + 2 points in B satisfying alternating condition and maximum absolute sum condition.
    :param B: Set of extreme_points
    :param d: The degree of the approximate polynomial
    :param r: The function for calculating the error
    """
    i = 0
    while i < len(B) - 1:
        if sign(r(B[i])) == sign(r(B[i + 1])):
            if abs(r(B[i])) < abs(r(B[i + 1])):
                B.pop(i)
            else:
                B.pop(i + 1)
        else:
            i += 1

    if len(B) > d + 3:
        # Calculate all |r(ti)| + |r(ti+1)| for i = 1, · · · , |B| − 1 and sort and store these values into the array T
        T = sorted([(abs(r(B[i])) + abs(r(B[i + 1])), i, i + 1) for i in range(len(B) - 1)])

    while len(B) > d + 2:
        if len(B) == d + 3:
            if abs(r(B[0])) < abs(r(B[-1])):
                B.pop(0)
            else:
                B.pop(-1)
        elif len(B) == d + 4:
            T.append((abs(r(B[0])) + abs(r(B[-1])), 0, len(B) - 1))
            T.sort()
            _, i1, i2 = T[0]
            B.pop(max(i1, i2))
            B.pop(min(i1, i2))
        else:
            if T[0][1] == 0 or T[0][2] == len(B) - 1:
                if T[0][1] == 0:
                    B.pop(0)
                else:
                    B.pop(-1)
            else:
                _, i1, i2 = T[0]
                B.pop(max(i1, i2))
                B.pop(min(i1, i2))

        # Recalculate T after removing elements
        T = sorted([(abs(r(B[i])) + abs(r(B[i + 1])), i, i + 1) for i in range(len(B) - 1)])

    return B


def second_deriv(coeffs):
    x_sym = sp.symbols('x', real=True)
    deriv_cheb = Chebyshev(coeffs)

    # Define diference function
    diff = lambda x: np.abs(x) / x

    # Define the difference function
    diff_func = deriv_cheb(x_sym) - diff(x_sym)
    # Calculate the second derivative of the difference function
    second_derivative = sp.diff(diff_func, x_sym, 2)
    return second_derivative


def first_deriv(coeffs):
    x_sym = sp.symbols('x', real=True)
    deriv_cheb = Chebyshev(coeffs)

    # Define diference function
    diff = lambda x: np.abs(x) / x

    # Define the difference function
    diff_func = deriv_cheb(x_sym) - diff(x_sym)
    # Calculate the second derivative of the difference function
    first_derivative = sp.diff(diff_func, x_sym)
    return first_derivative


def concavity(second_derivative, x_value):
    x_sym = sp.symbols('x', real=True)
    # Evaluate the second derivative at the point x_value
    # second_derivative_at_x = second_derivative.subs(x_sym, x_value)
    second_derivative_at_x = second_derivative(x_value)

    # Determine concavity
    if second_derivative_at_x < 0:
        return 1
    elif second_derivative_at_x > 0:
        return -1
    else:
        return 0


def filter_pos_min_neg_max(extremes, sec_derivative, diff):
    """
   Filter positive minimum and negative maximum points for a set of extreme points

   :param extremes: set of points to be filtered
   :param coeffs: coefficients of polynomial approximation p(x)
   :param diff: error function for p(x) and actual function

   :return: the polynomial coefficients, and an approximate maximum error associated with this approximation
   """

    filtered_extremes = []
    for extreme in extremes:
        concav = concavity(sec_derivative, extreme)
        if concav * sign(diff(extreme)) == 1:
            filtered_extremes.append(extreme)

    return filtered_extremes


def remez(func, n_degree: int, D: [(Decimal, Decimal)], max_iter: int = 10):
    """
    Approximate a function using the Remez algorithm

    :param func: a function (or lambda) f: X -> R
    :param n_degree: the degree of the polynomial to approximate the function f
    :param max_iter: Maximum number of iterations we try to converge for
    :param D: The list of intervals

    :return: the polynomial coefficients, and an approximate maximum error associated with this approximation
    """

    # initialize the node points
    n = n_degree + 1
    x_points = generate_points_in_ranges(D, n + 1)

    mean_error = float('inf')

    # for i in range(n + 1):
    #     A[i, n_degree + 1] = (-1) ** (i + 1)
    while True:
    # for _ in range(max_iter):
        A = mp.matrix(n + 1)
        coeffs = np.zeros(n_degree + 2)
        for i in range(n + 1):
            A[i, n_degree + 1] = (-1) ** (i + 1)
        # build the system
        vander = np.polynomial.chebyshev.chebvander(x_points, n_degree)

        for i in range(n + 1):
            for j in range(n):
                A[i, j] = vander[i, j]

        b = mp.matrix([mpf(func(x)) for x in x_points])

        # Convert mp.matrix A to numpy array
        A_np = np.array(A.tolist(), dtype=mpf)

        # Check if matrix A is singular or nearly singular
        try:
            det = mp.det(A_np)
        except TypeError:
            print("\'>=\' not supported between instances of \'NoneType\' and \'int\'")
            det = 0

        # cond_number = mp.cond(A_np)

        if det == 0:  # Adjust the threshold as needed
            print(f"Matrix A is singular or nearly singular. Determinant: {det}")
            print("Regularizing the matrix")
            # Regularize the matrix
            A_np_reg = A_np + np.eye(A_np.shape[0]) * 1e-20
            A = mp.matrix(A_np_reg)

        try:
            # l = mp.qr_solve (A, b)[0]
            l = mp.lu_solve(A, b)
        except ZeroDivisionError:
            print("LU Solve failed due to singular matrix.")
            # Handle or recover from the error
            quit(1)

        coeffs = l[:-1]

        coeffs_list = [coeffs[i, 0] for i in range(len(coeffs))]

        # debug = False
        #
        # num_samples = 3200
        #
        # if debug:
        #     r_i = lambda x: (np.polynomial.chebyshev.chebval(x, coeffs_list) - func(x))
        #     a = D[0][0]
        #     b = D[0][1]
        #     sample_points = np.linspace(a, b, num_samples)
        #     k = 0
        #     for i in range(num_samples - 1):
        #         x1, x2 = sample_points[i], sample_points[i + 1]
        #         if r_i(x1) * r_i(x2) < 0:
        #             k += 1
        #             print("gottem")
        #     print(k)


        # print("check_point_3")

        r_i_diff = lambda x: (np.polynomial.chebyshev.chebval(x, coeffs_list) - func(x))

        # This is a simplified version of the conditions for the extremes
        extremes = []

        temp_first_derivative = lambda x: mp.diff(r_i_diff, x)
        temp_second_derivative = lambda x: mp.diff(r_i_diff, x, 2)


        # second_derivative = second_deriv(coeffs)
        # first_derivative = first_deriv(coeffs)

        for lower, upper in D:
            extremes.append(mpf(lower))
            zero_derivatives = find_zeros(temp_first_derivative, lower, upper)
            filtered_zero_derivatives = filter_pos_min_neg_max(zero_derivatives, temp_second_derivative, r_i_diff)
            extremes.extend(filtered_zero_derivatives)

            extremes.append(mpf(upper))

        extremes.sort()

        # sets_of_extremes = find_combinations(extremes, n)
        # optimal_extremes = find_optimal_points(sets_of_extremes, r_i_diff)
        optimal_extremes = ChooseNewXpoints(extremes, n_degree, r_i_diff)
        errors = [abs(r_i_diff(i)) for i in optimal_extremes]
        mean_error = np.mean(errors)

        abs_errors = []
        # print("check_point_4")

        if min(errors) == 0:
          if np.max([abs(error - mean_error) for error in errors]) < approximation_param:
             break
        else:
            if (max(errors) - (min(errors))) / min(errors) < approximation_param:
                break

        # 
        x_points = optimal_extremes
        # return np.polynomial.chebyshev.cheb2poly(coeffs).tolist(), max(errors)

    return coeffs_list, max(errors)


# Please note that implementation is very specific for sign function
#  over two intervals
function = lambda x: mp.fabs(x) / x

# alpha = 20
# epsilon = mpf(2 ** (-20))

approximation_param = float(sys.argv[1])
D_string = sys.argv[2]

# Use regex to find all tuples
matches = re.findall(r'\(([^)]+)\)', D_string)

# Convert each match into a tuple of floats
D = [tuple(map(mpf, match.split(','))) for match in matches]
# print(D)

# degrees = [5, 5, 5, 7, 7, 7, 15]
# Degrees need to be odd for sign function
degrees_string = sys.argv[3][1:-1]
# print(degrees_string)
degrees_string_list = degrees_string.split(",")
degrees = [int(degree_string) for degree_string in degrees_string_list[:-1]] 

print("Approximating degrees: ", degrees)

i  = 1
for degree in degrees:
    print(f"Working on degree: {degree}")
    poly_coeffs, max_error = remez(function, degree, D, max_iter=10)
    D = [(mpf(-1), mpf(-1 + max_error)), (mpf(1 - max_error), mpf(1))]
    with open(f'../remez_outputs/coefficients_{degree}_{i}.txt', 'w') as f:
        for coeff in poly_coeffs:
            f.write(f"{coeff}\n")
    i+=1
    # print(D)

# D = [(mpf(-1.0), mpf(-0.99999999999999999993223145700552668175695437597355396)), (mpf(0.99999999999999999993223145700552668175695437597355396), mpf(1.0))]
# coeffs_list, max_error = remez(function, 5, D, max_iter=10)


# Plotting the polynomials
# x = np.linspace(-1, 1, 50)
# y_approx = np.polyval(poly_coeffs[::-1], x)
# y_exact = np.array([function(x_i) for x_i in x])
# plt.plot(x, y_exact, label='y exact')
# plt.plot(x, y_approx, label='x approximated')
# plt.legend()
# plt.title(f'$f(x)$ v. $P^*_{5}(x)$')
# plt.show()
