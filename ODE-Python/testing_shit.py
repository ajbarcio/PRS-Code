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

jac = np.array([[1,1]])
jacI = lin.pinv(jac)
print(jacI)
print(jacI.dot(jac))
print(jac.dot(jacI))

cmap = plt.get_cmap('hot')
print(cmap)

