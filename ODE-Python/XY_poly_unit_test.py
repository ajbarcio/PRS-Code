from stiffness_library import *
import numpy as np
import matplotlib.pyplot as plt

# use default spring for test
testSpring = Spring()

XCoeffs, YCoeffs = xy_poly(testSpring.pts, testSpring.XYArcLens)
smesh = np.linspace(0,testSpring.fullArcLength,101)
X = PPoly_Eval(smesh,XCoeffs)
Y = PPoly_Eval(smesh,YCoeffs)

dxds = PPoly_Eval(smesh, XCoeffs, deriv=1)
dyds = PPoly_Eval(smesh, YCoeffs, deriv=1)

plt.plot(X,Y)
plt.axis("equal")
plt.figure()
plt.plot(smesh, X)
plt.plot(smesh, dxds)
plt.figure()
plt.plot(smesh, Y)
plt.plot(smesh, dyds)
plt.show()

x = PPoly_Eval(0,XCoeffs)
print(x)