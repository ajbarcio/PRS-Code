import numpy as np
from scipy import optimize as opt
import scipy.integrate as intg
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
testCrsc1 = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = np.array([.0003, 0.000026, .0003]),
                       IcParamLens = np.array([0.5]))

# fuck this
testPath.get_crscRef(testCrsc)

# Initialize a spring and generate its shape (rootfinding method)
testSpring = Spring(testCrsc, materials.Maraging300Steel,
                    resolution=1000, name="20270730_spring")
testSpring.crsc.get_outer_geometry(testSpring.resl)

# Change this value to edit the index along the mesh at which we want to start
# the ODE method of thickness generation
smi = 0

# Set the initial value for the ODE
y0 = np.array([testSpring.crsc.la[smi], testSpring.crsc.lb[smi]])
print(y0)
# Integrate forward along the fixed mesh
geometry = fixed_rk4(testSpring.geo_ODE, y0, testSpring.ximesh[smi:])
print("about to try variable mesh")
python_res = intg.solve_ivp(testSpring.geo_ODE, (testSpring.ximesh[smi],
                                                 testSpring.ximesh[-1]), y0,
                                                 method='LSODA')
print(testSpring.singularityCounter)
print(testSpring.numericalSubCounter)
geometry_alt = python_res.y

plt.figure("Outer Profile Result Comparison")
plt.plot(testSpring.ximesh[smi:], testSpring.crsc.la[smi:], color='#cc1616', label="rootfinding la")
plt.plot(testSpring.ximesh[smi:], testSpring.crsc.lb[smi:], color='#3f0d80', label="rootfinding la")

plt.plot(testSpring.ximesh[smi:], geometry[0,:], "--", color='#16800d', label="ODE la")
plt.plot(testSpring.ximesh[smi:], geometry[1,:], "--", color='#c4c116', label="ODE lb")
plt.plot(python_res.t, geometry_alt[0,:], label="adaptive mesh ODE la")
plt.plot(python_res.t, geometry_alt[1,:], label="adaptive mesh ODE lb")


plt.ylim(0,.14)
plt.legend()

h = geometry[0,:]+geometry[1,:]
err = testSpring.crsc.h[smi:]-h

plt.figure("result error")
plt.plot(testSpring.ximesh[smi:], err, label="error")

# Empirically (numerically) find the derivatives of the thickness profile wrt s
dladxi_real = numerical_fixed_mesh_diff(testSpring.crsc.la, testSpring.ximesh)
dlbdxi_real = numerical_fixed_mesh_diff(testSpring.crsc.lb, testSpring.ximesh)

# Evaluate ODE to analytically find derivatives
dladxi_analytical = np.empty_like(dladxi_real)
dlbdxi_analytical = np.empty_like(dladxi_real)

i = 0
for value in testSpring.ximesh:
    dlalbdxi_analytical = testSpring.geo_ODE(value, [testSpring.crsc.la[i], testSpring.crsc.lb[i]])
    dladxi_analytical[i] = dlalbdxi_analytical[0]
    dlbdxi_analytical[i] = dlalbdxi_analytical[1]
    i+=1

plt.figure("derivative error")
plt.plot(testSpring.ximesh, dlbdxi_real-dlbdxi_analytical, label="dlbds error")
plt.plot(testSpring.ximesh, dladxi_real-dladxi_analytical, label="dlads error")
plt.ylim(-0.5,0.5)
plt.legend()

plt.figure("Derivative Result Comparison")
plt.plot(testSpring.ximesh, dladxi_real, label="real dlads")
plt.plot(testSpring.ximesh, dladxi_analytical, label="analytical dlads")
plt.plot(testSpring.ximesh, dlbdxi_real, label="real dlbds")
plt.plot(testSpring.ximesh, dlbdxi_analytical, label="analytical dlads")
plt.ylim(-.06,0.06)
plt.legend()

plt.figure("rn derivative accuracy check (does derivative act like a derivative)")
plt.plot(testSpring.ximesh, testSpring.path.get_rn(testSpring.ximesh), label="rn")
plt.plot(testSpring.ximesh, testSpring.path.get_drn(testSpring.ximesh), label="drn")
plt.plot(testSpring.ximesh, numerical_fixed_mesh_diff(testSpring.path.get_rn(testSpring.ximesh), testSpring.ximesh), label="drn numerical")
plt.ylim(-10,10)
plt.legend()

plt.show()