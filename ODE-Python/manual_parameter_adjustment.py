import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring import Spring


deg2rad = np.pi/180

################################################################################
# This script tests functions used to apply loading conditions to a spring geometry
################################################################################

# Fit the Gen 2 actuator form factor
R0 = 2.2/2
R3 = 5.9/2

# All these radii, angles, and Ic come from a previous version of this code and
# can be considered as arbitrary starting points

R1 = (R0+R3)/2+.26+.125
R2 = (R0+R3)/2-.25

fullAngle = 165

beta1 = fullAngle/3*deg2rad*.5

beta2 = 2*fullAngle/3*deg2rad*0.9
beta0 = fullAngle*deg2rad

# generate all the input arrays for adjustible parameters
# (treat material properties, fullParamLength guess, out of plane thickness,
#  and resolution as fixed)

inputRadii      = np.array([R0,R1,R2,R3])
inputBetaAngles = np.array([0,beta1,beta2,beta0])

Ics = np.array([0.008, 0.00025, 0.00025, 0.008])
IcLens=np.array([0.4, 0.66])

XYParamLens = np.array([0.333,0.667])

# Generate the spring:

curvedSpring = Spring(n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens)

testTorque = 4554.6
# select method
method = "smartGuess"
# for the smartGuess method, linear extrapolation for the BC forces determines
# an initial guess for the forces at full torque
if method=="smartGuess":
    print("using smart guess method")
    # get a ballpark initial guess
    SFGuess = curvedSpring.smart_initial_load_guess(testTorque,curvedSpring.deform_ODE)
    print("guess planar force vector:", SFGuess)
    # use that initial guess to deform the beam
    res, SF, divergeFlag, i = curvedSpring.wrapped_torque_deform(testTorque,curvedSpring.deform_ODE,SF=SFGuess)
    # compare the two guesses
    print("solution planar force vector:", SF)
    # print out results
    print("torque:", testTorque)
    print("angular deformation:", curvedSpring.dBeta/deg2rad)
    print("stiffness (lbf/deg)", testTorque/(curvedSpring.dBeta/deg2rad))
    print("max stress:", curvedSpring.maxStress)
    print("des stress:", curvedSpring.designStress)
    # if it didnt diverge consider the deform ation succesful
    deformBool = not(divergeFlag)
    # plot results
    curvedSpring.full_results(deformBool=deformBool, plotBool=True)
elif method=="slowRamp":
    print("using slow ramp method")
    # starting at an unloaded condition, calculate solutions for progressively increasing torque conditions
    res, SF, divergeFlag, i = curvedSpring.deform_by_torque_slowRamp(testTorque,curvedSpring.deform_ODE)
    deformBool = not(divergeFlag)
    print("solution planar force vector:", SF)
    print("angular deformation:", curvedSpring.dBeta/deg2rad)
    print("stiffness (lbf/deg)", testTorque/(curvedSpring.dBeta/deg2rad))
    # plot results
    curvedSpring.full_results(deformBool=deformBool, plotBool=True)

## output curves for use in solidworks
# add a "zero" column for z-values (needed by solidworks) and save as .txt files
A = np.hstack((curvedSpring.undeformedASurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedASurface))).T))
B = np.hstack((curvedSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedBSurface))).T))
np.savetxt("surfaces\\adjusted_A_surface.txt", A, delimiter=",", fmt='%f')
np.savetxt("surfaces\\adjusted_B_surface.txt", B, delimiter=",", fmt='%f')