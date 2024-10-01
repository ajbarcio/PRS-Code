import numpy as np
from scipy import optimize as opt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad
import matplotlib.pyplot as plt

from utils import numerical_fixed_mesh_diff

IR = 1
OR = 2

# Simplest outward spiral (center of curvature at origin)
testPath = PATHDEF.Empty(2,5,0,10)

testCrsc = CRSCDEF.Empty(testPath, 0.375)
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel, resolution=200)

err, res = testSpring.free_optimization(3000,SF=np.array([0,0,100]),
                                        RCs = np.array([0,testPath.x0,testPath.y0,.04,.05,10,90*deg2rad]))

# xorg
# yorg
rn = res[5,:]
la = res[3,:]
lb = res[4,:]
dgdxi = numerical_fixed_mesh_diff(res[0,:], testSpring.ximesh)
# no distortion
dgds  = dgdxi
stressInner = abs(testSpring.E*dgds*rn*(la/(rn-la)))
stressOuter = abs(testSpring.E*dgds*rn*(lb/(rn+lb)))

print(err)
plt.figure("deflection resutls")
plt.plot(res[0,:])
plt.figure("x, y results (deformed neutral surface)")
plt.plot(res[1,:],res[2,:])
plt.axis('equal')
plt.figure("la, lb results")
plt.plot(testSpring.ximesh, res[3,:])
plt.plot(testSpring.ximesh, res[4,:])
print(res[5,:])
plt.figure("rn results")
plt.plot(testSpring.ximesh,res[5,:])
plt.figure("alpha results")
plt.plot(testSpring.ximesh,res[6,:])
plt.figure("stress results")
plt.plot(testSpring.ximesh,stressInner)
plt.plot(testSpring.ximesh,stressOuter)



plt.show()