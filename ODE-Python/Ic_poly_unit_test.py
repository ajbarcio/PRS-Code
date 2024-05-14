from stiffness_library import *
import numpy as np

IcPts = np.array([.008,.001,.008])
IcArcLens = np.array([0,1,2])

IcCoeffs = Ic_poly(IcPts,IcArcLens)
print(IcCoeffs)
smesh = np.linspace(IcArcLens[0],IcArcLens[-1],101)
Ic = Ic_s(smesh,IcCoeffs)
dIcds = Ic_s(smesh,IcCoeffs,deriv=1)
d2Icds2 = Ic_s(smesh,IcCoeffs,deriv=2)
# d3Icds3 = Ic_s(smesh,IcCoeffs,deriv=3)
plt.plot(smesh,Ic)
plt.plot(smesh,dIcds)
# plt.plot(smesh,d2Icds2)
# plt.plot(smesh,d3Icds3)

plt.show()