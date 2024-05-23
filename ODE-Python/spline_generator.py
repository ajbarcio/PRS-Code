import numpy as np
from scipy.integrate import solve_ivp as variable_mesh_solve
import numpy.linalg as lin
import matplotlib.pyplot as plt
import time
from stiffness_library import *
from StatProfiler import SSProfile

ctrlLengths = [0,1,2,3]
pts = np.array([[0,0],[1,0],[2,0],[3,0]])
xCoeffs = xy_poly(pts, ctrlLengths)[0]
yCoeffs = xy_poly(pts, ctrlLengths)[1]
print(xCoeffs)
print(yCoeffs)
print(coord(0,xCoeffs), coord(0,yCoeffs))
print(coord(1,xCoeffs), coord(1,yCoeffs))
print(coord(2,xCoeffs), coord(2,yCoeffs))
print(coord(3,xCoeffs), coord(3,yCoeffs))

smesh = np.linspace(0,3,301)
xs = coord(smesh,xCoeffs)
ys = coord(smesh,yCoeffs)
plt.scatter(pts[:,0],pts[:,1])
plt.plot(xs,ys)
plt.show()