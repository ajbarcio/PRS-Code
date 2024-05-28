# from stiffness_library import *
import numpy as np
import matplotlib.pyplot as plt

from spring import Spring
from polynomials import xy_poly, PPoly_Eval

# use default spring for test
testSpring = Spring()

# generate polynomials in x and y
XCoeffs, YCoeffs = xy_poly(testSpring.pts, testSpring.XYParamLens)
# generate mesh over which to evaluate polynomial
mesh = np.linspace(0, testSpring.fullParamLength,101)
# generate x, y coordinates of surface
X = PPoly_Eval(mesh, XCoeffs)
Y = PPoly_Eval(mesh, YCoeffs)

# evaluate derivatives of x, y coordinates
dxdxi = PPoly_Eval(mesh, XCoeffs, deriv=1)
dydxi = PPoly_Eval(mesh, YCoeffs, deriv=1)

# plot generated neutral surface
plt.figure("Netural Surface")
plt.plot(X,Y)
plt.axis("equal")
# plot x polynomial with derivative to ensure it makes sense
plt.figure("X Poly")
plt.plot(mesh, X, label = "X")
plt.plot(mesh, dxdxi, label = "dxdxi")
plt.legend()
# plot y polynomial with derivative to ensure it makes sense
plt.figure("Y Poly")
plt.plot(mesh, Y, label = "Y")
plt.plot(mesh, dydxi, label = "dydxi")
plt.legend()
plt.show()
