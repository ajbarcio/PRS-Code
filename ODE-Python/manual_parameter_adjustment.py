import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring import Spring
from materials import Maraging300Steel

deg2rad = np.pi/180

################################################################################
# This script tests functions used to apply loading conditions to a spring geometry
################################################################################

# Fit the Gen 2 actuator form factor
R0 = 2.2/2
R3 = 5.9/2

# All these radii, angles, and Ic come from a previous version of this code and
# can be considered as arbitrary starting points

R1 = (R0+R3)/2+.26+.125+.125
R2 = (R0+R3)/2-.25+.125

fullAngle = 163

beta1 = fullAngle/3*deg2rad*.5*1.3+15*deg2rad

beta2 = 2*fullAngle/3*deg2rad*0.9*.9+30*deg2rad
beta0 = fullAngle*deg2rad

# generate all the input arrays for adjustible parameters
# (treat material properties, fullParamLength guess, out of plane thickness,
#  and resolution as fixed)

inputRadii      = np.array([R0,R1,R2,R3])
inputBetaAngles = np.array([0,beta1,beta2,beta0])

Ics = np.array([0.008*.85, 0.00025*.25, 0.00025*.25, 0.008*1.15])
IcLens=np.array([0.4*1.2, 0.667*0.8])

XYParamLens = np.array([0.333,0.667*1.02])

# Generate the spring:

curvedSpring = Spring(Maraging300Steel(), n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens, XYParamLens=XYParamLens,
                             name="manual_spring")

testTorque = 4554.6
# select method
method = "smartGuess"
# for the smartGuess method, linear extrapolation for the BC forces determines
# an initial guess for the forces at full torque
if method=="smartGuess":
    print("using smart guess method")
    # use smart guess method to solve spring
    res, SF, divergeFlag, i = curvedSpring.deform_by_torque_smartGuess \
                             (testTorque,curvedSpring.deform_ODE,printBool=True)
    # if it didnt diverge consider the deform ation succesful
    deformBool = not(divergeFlag)
    # plot results
    curvedSpring.full_results(deformBool=deformBool, plotBool=True)
    # print out results
    print("solution planar force vector:", SF)
    print("torque:", testTorque)
    print("angular deformation:", curvedSpring.dBeta/deg2rad)
    print("stiffness (lbf/deg)", testTorque/(curvedSpring.dBeta/deg2rad))
    print("ideal stiffness:   ", testTorque/(5))
    print("max stress:", curvedSpring.maxStress)
    print("des stress:", curvedSpring.designStress)
elif method=="slowRamp":
    print("using slow ramp method")
    # starting at an unloaded condition, calculate solutions for progressively increasing torque conditions
    res, SF, divergeFlag, i = curvedSpring.deform_by_torque_slowRamp(testTorque,curvedSpring.deform_ODE,torqueResolution=50)
    deformBool = not(divergeFlag)
    print("solution planar force vector:", SF)
    print("angular deformation:", curvedSpring.dBeta/deg2rad)
    print("stiffness (lbf/deg)", testTorque/(curvedSpring.dBeta/deg2rad))
    # plot results
    curvedSpring.full_results(deformBool=deformBool, plotBool=True)

curvedSpring.full_results(deformBool=False)
## output curves for use in solidworks
curvedSpring.surfaces_to_txt()

plt.show()