import numpy as np
import numpy.linalg as lin

from matplotlib import pyplot as plt

from CRSCDEF import Constant_Ic
from PATHDEF import LinearRnSpiral, RadiallyEndedPolynomial

# path = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3))
path = LinearRnSpiral(2, 2*np.pi/3, initialRadius=1, finalRadius=5)
# path = RadiallyEndedPolynomial(2, 6)
crsc = Constant_Ic(path, 0.1, Ic0=0.125)

ximesh = np.linspace(0,path.arcLen,101)
targetIcs = crsc.get_Ic(ximesh)
lalb      = crsc.get_lalb(ximesh)
la = lalb[0,:]
lb = lalb[1,:]
print(la[-5:])
print(lb[-5:])
rn = path.get_rn(ximesh)
actualIcs = crsc.t*rn/2*(lb**2-la**2)

plt.plot(ximesh, targetIcs)
plt.plot(ximesh, actualIcs,"--")
plt.show()