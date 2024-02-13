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

jac = np.array([[1,2],[2,1]])
jacInv = lin.inv(jac)
print(jacInv)
vector = np.array([1,1])
print(jacInv.dot(vector))
print(vector.dot(jacInv))
