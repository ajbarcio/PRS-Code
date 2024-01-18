import time
import math
import numpy as np
import scipy as sp
from scipy import integrate
from scipy import optimize as op
import numpy.linalg as lin
import scipy.linalg as slin
import matplotlib.pyplot as plt
import matplotlib.transforms as tfm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *
from openpyxl import *


p0 = np.array([0,0])
p1 = np.array([0,1])
p2 = np.array([0,2])
p3 = np.array([0,3])

pts = np.array([p0,p1,p2,p3])

u = np.linspace(0,1,101)
u2 = u - np.floor(u)
xyc = np.empty([len(u), 2])

print(u)

for i in range(len(u)):
        t = u2[i]
        U = np.array([t**3, t**2, t, 1])
        D = np.array([[-1,3,-3,1],[3,-6,3,0],[-3,3,0,0],[1,0,0,0]])
        P = np.array([pts[0],pts[1],pts[2],pts[3]])
        xyc[i] = U.dot(D).dot(P)
plt.plot(u,xyc[:,1])
plt.show()