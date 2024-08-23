import numpy as np
from scipy import optimize as opt
from matplotlib import pyplot as plt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

from utils import fixed_rk4, numerical_fixed_mesh_diff

# Standard procedure: define a path, then define a thickness profile using the
#                     length of that path, then pass both of those objects into
#                     the overall spring object that contains deformation
#                     methods

IR = 1.3
OR = 2.013

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
                                       radii = np.array([IR,(IR+OR)/2*1.15,OR]),
                                       ffradii = np.array([1.14, 2.113]),
                                       alphaAngles = np.array([45,45])*deg2rad,
                                       betaAngles = np.array([0,87.5,175])*deg2rad,
                                       XYFactors = np.array([0.5]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = np.array([.002, .000026, .000026, .0003]),
                       IcParamLens = np.array([.50, .6]))
testCrsc2 = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = np.array([.0003, 0.000026, .0003]),
                       IcParamLens = np.array([0.5]))

# fuck this
testPath.get_crscRef(testCrsc2)


testSpring = Spring(testCrsc2, materials.Maraging300Steel, name="20270730_spring")
testSpring.crsc.get_outer_geometry(testSpring.resl)

smi = 0

y0 = np.array([testSpring.crsc.la[smi], testSpring.crsc.lb[smi]])

geometry = fixed_rk4(testSpring.geo_ODE, y0, testSpring.ximesh[smi:])
# print(geometry)

plt.figure()
plt.plot(testSpring.ximesh[smi:], geometry[0,:])
plt.plot(testSpring.ximesh[smi:], geometry[1,:])
plt.plot(testSpring.ximesh[smi:], testSpring.crsc.la[smi:])
plt.plot(testSpring.ximesh[smi:], testSpring.crsc.lb[smi:])

dladxi_real = numerical_fixed_mesh_diff(testSpring.crsc.la, testSpring.ximesh)
dlbdxi_real = numerical_fixed_mesh_diff(testSpring.crsc.lb, testSpring.ximesh)

dladxi_analytical = np.empty_like(dladxi_real)
dlbdxi_analytical = np.empty_like(dladxi_real)

i = 0
for value in testSpring.ximesh:
    dlalbdxi_analytical = testSpring.geo_ODE(value, [testSpring.crsc.la[i], testSpring.crsc.lb[i]])
    dladxi_analytical[i] = dlalbdxi_analytical[0]
    dlbdxi_analytical[i] = dlalbdxi_analytical[1]
    i+=1

# print(dladxi_analytical)
# print(dlbdxi_analytical)
plt.figure()
plt.plot(testSpring.ximesh, dlbdxi_real-dlbdxi_analytical, label="b error")
plt.plot(testSpring.ximesh, dladxi_real-dladxi_analytical, label="a error")
plt.ylim(-0.5,0.5)
plt.legend()
plt.figure()
plt.plot(testSpring.ximesh, dladxi_real, label="real a")
plt.plot(testSpring.ximesh, dladxi_analytical, label="analytical a")
plt.plot(testSpring.ximesh, dlbdxi_real, label="real b")
plt.plot(testSpring.ximesh, dlbdxi_analytical, label="analytical b")

allRn = testSpring.path.get_rn(testSpring.ximesh)

plt.plot(testSpring.ximesh, allRn, label="rn")
plt.ylim(-.06,0.06)
plt.legend()

plt.figure()
plt.plot(testSpring.ximesh, testSpring.path.get_rn(testSpring.ximesh), label="rn")
plt.plot(testSpring.ximesh, testSpring.path.get_drn(testSpring.ximesh), label="drn")
plt.plot(testSpring.ximesh, numerical_fixed_mesh_diff(testSpring.path.get_rn(testSpring.ximesh), testSpring.ximesh), label="drn numerical")
plt.ylim(-10,10)
plt.legend()

plt.show()