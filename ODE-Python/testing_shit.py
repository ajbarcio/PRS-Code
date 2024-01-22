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

u = np.linspace(0,3,3001)
gammaSol = np.zeros((u.shape[0], 2))
fuck = np.array([1, 1])
gammaSol[-1] = fuck
print(gammaSol)