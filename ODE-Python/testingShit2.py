from stiffness_library import *
import numpy as np

IcPts = np.array([3,1,3])
IcArcLens = np.array([0,1,2])

IcCoeffs = Ic_poly(IcPts,IcArcLens)
print(IcCoeffs)
smesh = np.linspace(IcArcLens[0],IcArcLens[-1],101)
Ic = PPoly_Eval(smesh,IcCoeffs)
plt.plot(smesh,Ic)
plt.show()