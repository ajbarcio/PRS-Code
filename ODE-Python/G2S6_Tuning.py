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

IR = springData.loc['Size 6','IR lim (in)']
OR = springData.loc['Size 6','OR lim (in)']
testTorque = springData.loc['Size 6','Max Torque (in.lbs)']
# IR = 1.3
# OR = 2.013

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
                                       radii = np.array([IR,(IR+OR)/2*1.15,OR]),
                                       # ffradii = np.array([1.14, 2.113]),
                                       ffradii = np.array([1.14, 2.5]),
                                       # alphaAngles = np.array([45,45])*deg2rad,
                                       alphaAngles = np.array([55,0])*deg2rad,
                                    #    betaAngles = np.array([0,87.5,175])*deg2rad,
                                       betaAngles = np.array([0,95,175])*deg2rad,
                                       XYFactors = np.array([0.5]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                                        # IcPts = np.array([.002, .000026, .000026, .0003]),
                       IcPts = np.array([.0038, .000025, .000015, .00001, .00015])*2,
                       IcParamLens = np.array([.45, .50, .65]))
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Titanium5, name="G2S6_spring_Titanium")

# Actually try to deform it
# res, SF, divergeFlag, i = testSpring.deform_by_torque(-4549,
#                                                       testSpring.deform_ODE,
#                                                       SF=np.array([0,0,-4549]))
# testSpring.plot_deform(showBool=False)
res, SF, divergeFlag, i = testSpring.deform_by_torque(testTorque,
                                                      testSpring.deform_ODE,
                                                      SF=np.array([0,0,4549]))

print(testSpring.solnerr)
print(SF)
testSpring.plot_deform()

testSpring.export_surfaces()
