import numpy as np
from numpy import linalg as lin
from scipy import optimize as opt
from matplotlib import pyplot as plt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

from utils import fixed_rk4, numerical_fixed_mesh_diff

import sympy as sp

x = np.array([[0],[float('nan')]])
print(x)
print(np.isfinite(x))
print(np.isfinite(x).any())
print(np.invert(np.isfinite(x)))
print(np.invert(np.isfinite(x)).any())
x = np.array([[0],[5]])
print(x)
print(np.isfinite(x))
print(np.isfinite(x).any())
print(np.invert(np.isfinite(x)))
print(np.invert(np.isfinite(x)).any())


# x = sp.symbols('x')

# func = (x**2-4)/(3*x+6)

# print(sp.limit(func, x, -2))

# coeffs = np.array([1,2.2,3,4,5.2])
# polynomial = sp.Poly(coeffs, x).as_expr()
# differential = sp.diff(polynomial, x)
# func = sp.cos(x)
# print(sp.diff(polynomial, x))
# print(differential)

# # polynomial = polynomial.as_expr()
# # differential = differential.as_expr()

# input = np.array([1,2,3,4])

# polyLambda = sp.lambdify(x, polynomial, "numpy")
# diffLambda = sp.lambdify(x, differential)
# funcLambda = sp.lambdify(x, func,       "numpy")

# print(polyLambda(input))
# print(diffLambda(input))
# print(funcLambda(input))