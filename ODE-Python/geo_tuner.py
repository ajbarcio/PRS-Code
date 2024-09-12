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
gainValues = np.linspace(0,2,201)
maxErr = np.empty(len(gainValues))
# Integrate forward along the fixed mesh
i = 0
for gain in gainValues:
    print(gain)
    testSpring.geoFeedbackGain = -gain
    geometry = fixed_rk4(testSpring.geo_ODE, y0, testSpring.ximesh[smi:])

    h = geometry[0,:]+geometry[1,:]
    err = abs(testSpring.crsc.h[smi:]-h)
    errMax = np.max(err)
    maxErr[i] = errMax
    i+=1
plt.plot(gainValues, maxErr)
plt.show()