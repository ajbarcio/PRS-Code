import numpy as np
from numpy import linalg as lin
from scipy import optimize as opt
from matplotlib import pyplot as plt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

from utils import fixed_rk4, numerical_fixed_mesh_diff


a = 1
r = float('inf')
I = 1/12

A = np.array([[-a,a],[1/(r-a)-1/r,1/(r+a)+1/r]])
print(A)
print(lin.pinv(A))