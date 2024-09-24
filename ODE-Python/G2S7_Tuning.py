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

IR = springData.loc['Size 7','IR lim (in)']
OR = springData.loc['Size 7','OR lim (in)']
testTorque = springData.loc['Size 7','Max Torque (in.lbs)']

preDefT           = .375
preDefIc          = np.array([.008, .0001, .00005, .003])
preDefAlphaAngles = np.array([0,0])*deg2rad
preDefBetaRange   = 168

rootThk  = np.cbrt(12*preDefIc[0]/preDefT)
tipThk   = np.cbrt(12*preDefIc[-1]/preDefT)
thks = np.array(rootThk,tipThk)
offsets = thks*np.sin(preDefAlphaAngles)/2
offsets[-1] = -offsets[-1]

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
            radii = np.array([IR+offsets[0],(IR+OR)/2*1.1,(IR+OR)/2*1.1,OR+offsets[1]]),
            ffradii = np.array([IR, OR]),
            alphaAngles = preDefAlphaAngles,
            betaAngles = np.array([0,
                                    preDefBetaRange/3*.9,
                                    2*preDefBetaRange/3*1.1,
                                    preDefBetaRange])*deg2rad,
            XYFactors = np.array([0.3333,0.6666]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = preDefIc,
                       IcParamLens = np.array([0.4,0.6]),
                       outPlaneThickness = preDefT)
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel, name="G2S7_spring")

# Actually try to deform it
# res, SF, divergeFlag, i = testSpring.deform_by_torque(-4549,
#                                                       testSpring.deform_ODE,
#                                                       SF=np.array([0,0,-4549]))
# testSpring.plot_deform(showBool=False)
res, SF, divergeFlag, i = testSpring.deform_by_torque(testTorque,
                                                      testSpring.deform_ODE,
                                                      SF=np.array([0,0,testTorque]))

testSpring.linearity_plot(testTorque,200)

print(testSpring.solnerr)
print(SF)
testSpring.plot_deform()

testSpring.export_surfaces()
