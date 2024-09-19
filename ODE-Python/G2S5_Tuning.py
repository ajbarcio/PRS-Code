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

# Some inter-effecting parameters get defined here, ahead of time
preDefT           = .3
preDefIc          = np.array([.0008, .000012, .000012, .00085])
preDefAlphaAngles = np.array([75,30])*deg2rad
preDefBetaRange   = 175

rootThk  = np.cbrt(12*preDefIc[0]/preDefT)
tipThk   = np.cbrt(12*preDefIc[-1]/preDefT)
thks = np.array(rootThk,tipThk)
offsets = thks*np.sin(preDefAlphaAngles)/2
offsets[-1] = -offsets[-1]

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
            radii = np.array([IR+offsets[0],(IR+OR)/2*1.05,(IR+OR)/2*0.9,OR+offsets[1]]),
            ffradii = np.array([IR, OR]),
            alphaAngles = preDefAlphaAngles,
            betaAngles = np.array([0,preDefBetaRange/3*.85,
                                    2*preDefBetaRange/3,
                                    preDefBetaRange])*deg2rad,
            XYFactors = np.array([0.333,0.666]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = preDefIc,
                       IcParamLens = np.array([0.51, .7]),
                       outPlaneThickness = preDefT)
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel, name="G2S5_spring")

# Deform in a (hopefully) robust way
res, SF, divergeFlag, i = testSpring.deform_by_torque_predict_forces(testTorque,
                                                      testSpring.deform_ODE,
                                                      breakBool=True)
if divergeFlag:
    testSpring.detailed_deform_regression(testTorque, testSpring.deform_ODE,resl=10,degree=1)

print("total deformation:", testSpring.dBeta/deg2rad, u'\N{DEGREE SIGN}')
print(testSpring.solnerrVec)
print(testSpring.solnerr)
# print(SF)

print("Trying to optimize for deflection")

testSpring.plot_deform()

testSpring.export_surfaces()
