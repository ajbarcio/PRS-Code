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
testTorque = springData.loc['Size 6','Max Torque (in.lbs)']*0.5
# IR = 1.3
# OR = 2.013

print(IR, OR)

pren=2
preDefT           = .375
preDefIc          = np.array([.0028, .00002, .00002, .00015])
preDefAlphaAngles = np.array([35,25])*deg2rad
preDefBetaRange   = 175

rootThk  = np.cbrt(12*preDefIc[0]/preDefT)
tipThk   = np.cbrt(12*preDefIc[-1]/preDefT)
thks = np.array(rootThk,tipThk)
offsets = thks*np.sin(preDefAlphaAngles)/2
offsets[-1] = -offsets[-1]

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=pren, fullParamLength=3,
            radii = np.array([IR+offsets[0],(IR+OR)/2,(IR+OR)/2*.85,OR+offsets[1]]),
            ffradii = np.array([IR, OR]),
            alphaAngles = preDefAlphaAngles,
            betaAngles = np.array([0,preDefBetaRange*.3,preDefBetaRange*.75,preDefBetaRange])*deg2rad,
            XYFactors = np.array([0.3,0.6]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                                        IcPts = preDefIc,
                                        IcParamLens = np.array([.45, .5]),
                                        outPlaneThickness = preDefT)

# testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
#                                        radii = np.array([IR,(IR+OR)/2*1.15,OR]),
#                                        # ffradii = np.array([1.14, 2.113]),
#                                        ffradii = np.array([1.14, 2.5]),
#                                        # alphaAngles = np.array([45,45])*deg2rad,
#                                        alphaAngles = np.array([55,0])*deg2rad,
#                                     #    betaAngles = np.array([0,87.5,175])*deg2rad,
#                                        betaAngles = np.array([0,95,175])*deg2rad,
#                                        XYFactors = np.array([0.5]))
# testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
#                                         # IcPts = np.array([.002, .000026, .000026, .0003]),
#                        IcPts = np.array([.0038, .000025, .000015, .00001, .00015])*2,
#                        IcParamLens = np.array([.45, .50, .65]))
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Titanium5, name="G2S6_spring_Titanium")

# Deform in a (hopefully) robust way
res, SF, divergeFlag, i = testSpring.deform_by_torque_predict_forces(testTorque,
                                                      testSpring.deform_ODE,
                                                      breakBool=True)
if divergeFlag:
    testSpring.detailed_deform_regression(testTorque, testSpring.deform_ODE,resl=10,degree=1)

print("total deformation:", testSpring.dBeta/deg2rad, u'\N{DEGREE SIGN}')
print(testSpring.solnerrVec)
print(testSpring.solnerr)
print(SF)

testSpring.plot_deform()

testSpring.export_surfaces()
