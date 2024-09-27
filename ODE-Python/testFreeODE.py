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


testPath = PATHDEF.Empty(2,1,0,-2,-2)

testCrsc = CRSCDEF.Empty(testPath, 0.375)
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel, resolution=1000)

err, res = testSpring.free_optimization(3000,SF=np.array([0,0,3000]),
                                        RCs = np.array([0,0,10,testPath.x0,testPath.y0,0.04,0.05]))

print(err)
plt.plot(res[1,:],res[2,:])
plt.axis('equal')
plt.show()