import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring_old import Spring
from materials import Maraging300Steel

deg2rad = np.pi/180

################################################################################
# This script tests functions used to apply loading conditions to a spring geometry
################################################################################

# Fit the Gen 2 actuator form factor
R0 = 2.0/2
R3 = 5.0/2

# All these radii, angles, and Ic come from a previous version of this code and
# can be considered as arbitrary starting points

R1 = (R0+R3)/2+.26+.125+.125
R2 = (R0+R3)/2-.25+.125

fullAngle = 165

beta1 = fullAngle/3*deg2rad*.5*1.2

beta2 = 2*fullAngle/3*deg2rad*0.9*.9
beta0 = fullAngle*deg2rad

# generate all the input arrays for adjustible parameters
# (treat material properties, fullParamLength guess, out of plane thickness,
#  and resolution as fixed)

inputRadii      = np.array([R0,R1,R2,R3])
inputBetaAngles = np.array([0,beta1,beta2,beta0])

Ics = np.array([0.008*.85, 0.00025*.85, 0.00025*1.15, 0.008*1.15])
IcLens=np.array([0.4*1.2, 0.667*0.8])

XYParamLens = np.array([0.333,0.667*1.02])

# Generate the spring:

curvedSpring = Spring(Maraging300Steel(), n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens, XYParamLens=XYParamLens,
                             name="manual_spring")
R1 = R3/3
R2 = R3/3*2
inputRadii      = np.array([R0,R1,R2,R3])
straightThickness = 0.25
outThick = 0.375
IcStraight = outThick*straightThickness**3/12
# print(IcStraight)

straightSpring = Spring(Maraging300Steel(), n=6, radii=inputRadii,
                                        betaAngles=np.array([0,0,0,0]),
                                        outPlaneThickness=outThick,
                                        IcPts=np.array([IcStraight,IcStraight]),
                                        IcParamLens = np.array([0,1]),
                                        XYParamLens=np.array([1/3.0,2/3.0]),
                                        name="testStraightSpring")

# testTorque = 4554.6
# # select method
# method = "smartGuess"
# # for the smartGuess method, linear extrapolation for the BC forces determines
# # an initial guess for the forces at full torque
# if method=="smartGuess":
#     print("using smart guess method")
#     # get a ballpark initial guess
#     SFGuess = curvedSpring.smart_initial_load_guess(testTorque,curvedSpring.deform_ODE)
#     print("guess planar force vector:", SFGuess)
#     # use that initial guess to deform the beam
#     res, SF, divergeFlag, i = curvedSpring.wrapped_torque_deform(testTorque,curvedSpring.deform_ODE,SF=SFGuess)

#     # if it didnt diverge consider the deform ation succesful
#     deformBool = not(divergeFlag)
#     # plot results
#     curvedSpring.full_results(deformBool=deformBool, plotBool=True)

#     # compare the two guesses
#     print("solution planar force vector:", SF)
#     # print out results
#     print("torque:", testTorque)
#     print("angular deformation:", curvedSpring.dBeta/deg2rad)
#     print("stiffness (lbf/deg)", testTorque/(curvedSpring.dBeta/deg2rad))
#     print("ideal stiffness:   ", testTorque/(5))
#     print("max stress:", curvedSpring.maxStress)
#     print("des stress:", curvedSpring.designStress)
# elif method=="slowRamp":
#     print("using slow ramp method")
#     # starting at an unloaded condition, calculate solutions for progressively increasing torque conditions
#     res, SF, divergeFlag, i = curvedSpring.deform_by_torque_slowRamp(testTorque,curvedSpring.deform_ODE,torqueResolution=50)
#     deformBool = not(divergeFlag)
#     print("solution planar force vector:", SF)
#     print("angular deformation:", curvedSpring.dBeta/deg2rad)
#     print("stiffness (lbf/deg)", testTorque/(curvedSpring.dBeta/deg2rad))
#     # plot results
#     curvedSpring.full_results(deformBool=deformBool, plotBool=True)

# curvedSpring.full_results(deformBool=False)
# ## output curves for use in solidworks
# # add a "zero" column for z-values (needed by solidworks) and save as .txt files
# A = np.hstack((curvedSpring.undeformedASurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedASurface))).T))
# B = np.hstack((curvedSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedBSurface))).T))
# np.savetxt("surfaces\\adjusted_A_surface.txt", A, delimiter=",", fmt='%f')
# np.savetxt("surfaces\\adjusted_B_surface.txt", B, delimiter=",", fmt='%f')

# plt.show()


testTorque = 4554.6
# select method
method = "smartGuess"
# for the smartGuess method, linear extrapolation for the BC forces determines
# an initial guess for the forces at full torque
if method=="smartGuess":
    print("using smart guess method")
    # get a ballpark initial guess
    SFGuess = straightSpring.smart_initial_load_guess(testTorque,straightSpring.deform_ODE)
    print("guess planar force vector:", SFGuess)
    # use that initial guess to deform the beam
    res, SF, divergeFlag, i = straightSpring.wrapped_torque_deform(testTorque,straightSpring.deform_ODE,SF=SFGuess)

    # if it didnt diverge consider the deform ation succesful
    deformBool = not(divergeFlag)
    # plot results
    straightSpring.full_results(deformBool=deformBool, plotBool=True)

    # compare the two guesses
    print("solution planar force vector:", SF)
    # print out results
    print("torque:", testTorque)
    print("angular deformation:", straightSpring.dBeta/deg2rad)
    print("stiffness (lbf/deg)", testTorque/(straightSpring.dBeta/deg2rad))
    print("ideal stiffness:   ", testTorque/(5))
    print("max stress:", straightSpring.maxStress)
    print("des stress:", straightSpring.designStress)
elif method=="slowRamp":
    print("using slow ramp method")
    # starting at an unloaded condition, calculate solutions for progressively increasing torque conditions
    res, SF, divergeFlag, i = straightSpring.deform_by_torque_slowRamp(testTorque,straightSpring.deform_ODE,torqueResolution=50)
    deformBool = not(divergeFlag)
    print("solution planar force vector:", SF)
    print("angular deformation:", straightSpring.dBeta/deg2rad)
    print("stiffness (lbf/deg)", testTorque/(straightSpring.dBeta/deg2rad))
    # plot results
    straightSpring.full_results(deformBool=deformBool, plotBool=True)

straightSpring.full_results(deformBool=False)
## output curves for use in solidworks
# add a "zero" column for z-values (needed by solidworks) and save as .txt files
A = np.hstack((straightSpring.undeformedASurface,np.atleast_2d(np.zeros(len(straightSpring.undeformedASurface))).T))
B = np.hstack((straightSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(straightSpring.undeformedBSurface))).T))
np.savetxt("surfaces\\adjusted_A_surface.txt", A, delimiter=",", fmt='%f')
np.savetxt("surfaces\\adjusted_B_surface.txt", B, delimiter=",", fmt='%f')

plt.show()