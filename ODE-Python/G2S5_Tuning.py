import numpy as np
from scipy import optimize as opt

import pandas as pd

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

# Standard procedure: define a path, then define a thickness profile using the
#                     length of that path, then pass both of those objects into
#                     the overall spring object that contains deformation
#                     methods

springData = pd.read_excel('Spring_Constraints.ods', engine='odf', index_col=0)

IR = springData.loc['Size 5','IR lim (in)']
OR = springData.loc['Size 5','OR lim (in)']
testTorque = springData.loc['Size 5','Max Torque (in.lbs)']

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
                                       radii = np.array([IR,(IR+OR)/2*1.15,OR]),
                                       ffradii = np.array([IR, OR]),
                                       alphaAngles = np.array([0,5])*deg2rad,
                                       betaAngles = np.array([0,87.5,175])*deg2rad,
                                       XYFactors = np.array([0.5]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = np.array([.0003, .000026, .000026, .0005])*2,
                       IcParamLens = np.array([.4, .6]))
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel, name="G2S5_spring")

# Deform in a (hopefully) robust way
res, SF, divergeFlag, i = testSpring.deform_by_torque_predict_forces(-testTorque,
                                                      testSpring.deform_ODE,
                                                      breakBool=False)

print(testSpring.solnerrVec)
print(testSpring.solnerr)
# print(SF)
testSpring.plot_deform()

testSpring.export_surfaces()
