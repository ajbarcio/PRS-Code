# from stiffness_library import *
import numpy as np
import matplotlib.pyplot as plt
import numpy.linalg as lin

from spring import Spring
from polynomials import xy_poly, PPoly_Eval

# TODO: Add condition check to this test
#       (Sort of done, checks that zeroth-order constraints are kept)

###############################################################################
# This script tests the functions for creating x-y polynomials for the neutral
#      surface
###############################################################################

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

for i in range(len(testSpring.XYParamLens)):
    diffX = PPoly_Eval(testSpring.XYParamLens[i],XCoeffs)-testSpring.pts[i,0]
    diffY = PPoly_Eval(testSpring.XYParamLens[i],YCoeffs)-testSpring.pts[i,1]
    diff = lin.norm([diffX,diffY])
    assert(np.isclose(diff,0))
    print(i, diff)
print("above quantities should all be near or equal to zero")


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
