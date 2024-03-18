import numpy as np
import scipy as sp
from scipy import optimize as op
from scipy import interpolate as spint
import numpy.linalg as lin
import scipy.linalg as slin
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *
from openpyxl import *
import random

# jac = np.array([[1,2],[2,1]])
# jacInv = lin.inv(jac)
# print(jacInv)
# vector = np.array([1,1])
# print(jacInv.dot(vector))
# print(vector.dot(jacInv))

# s = np.array([1,2,3,4])
# t = np.array([4, 3, 2, 1])
# # U = np.array([s**4, s**3, s**2, s, 1])

# print(s.dot(np.identity(4)))
# print(s)
# print(np.identity(4).dot(s))
# print(np.identity(4)*(s))

def function(a, b):
    return a, b, 3

print(function(1,2)[0])

var = [0, 1]
for i in range(7):
    tempVar = var
    tempVar[0] = tempVar[0]+i
    print("tempVar:", tempVar)
    print("var:",var)