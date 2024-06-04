import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt\

from spring import Spring
from materials import * # This module ONLY includes materials objects so I think this is safe

deg2rad = np.pi/180

################################################################################
# This script tests functions used to apply loading conditions to a
# spring geometry
################################################################################

# Fit the Gen 2 actuator form factor
R0 = 2.2/2
R3 = 5.9/2

# All these radii, angles, and Ic come from a previous version of this code and
# can be considered as arbitrary starting points

R1 = (R0+R3)/2+.26
R2 = (R0+R3)/2-.25

fullAngle = 145

beta1 = fullAngle/3*deg2rad*.5
beta2 = 2*fullAngle/3*deg2rad
beta0 = fullAngle*deg2rad

Ics = np.array([0.008, 0.000025, 0.008])

# Generate three springs for test and comparison:
# curvedSpring:        An arbitrarily notional idea of what a spring should
#                         look like
# straightSpring:      A perfectly straight beam with constant cross-section
# quasiStraightSpring: A mostly-straight beam that should behave similarly to a
#                         straight beam
# THE ONLY DIFFERENCE BETWEEN straight AND quasiStraight IS THE BETA ANGLES

allSpringsMaterial = Maraging300Steel()
print(allSpringsMaterial.yieldStress)
curvedSpring = Spring(allSpringsMaterial, n = 2, radii=np.array([R0,R1,R2,R3]),
                             betaAngles=np.array([0,beta1,beta2,beta0]),
                             IcPts=Ics,
                             IcParamLens=np.array([0.6]), resolution=500)
straightSpring = Spring(allSpringsMaterial, n = 1, fullParamLength = 6, radii = np.array([1,3,5,7]),
                        betaAngles=np.array([0,0,0,0])*deg2rad,
                        IcPts = np.array([0.03125, 0.03125, 0.03125]), resolution = 200)
quasiStraightSpring = Spring(allSpringsMaterial, n = 1, fullParamLength = 6, radii = np.array([1,3,5,7]),
                        betaAngles=np.array([0,1,2,3])*deg2rad,
                        IcPts = np.array([0.03125, 0.03125, 0.03125]), resolution = 200)

############################ Unit Test 0 ############################
# Make sure that for the degenerate case, (beams with 0 torque)        #
# the expected values are returned (no deflection at all)              #
########################################################################

print("#################### UNIT TEST 0: Checking 0-torque case ####################")
# perform a single pass of solving the deformation ODE under 0 load for each spring
err, res = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([0, 0, 0]), 0)
# if the spring does not deform, error should be near 0
assert(np.isclose(lin.norm(err),0))

err, res = quasiStraightSpring.forward_integration(quasiStraightSpring.deform_ODE, np.array([0, 0, 0]), 0)
assert(np.isclose(lin.norm(err),0))

err, res = curvedSpring.forward_integration(curvedSpring.deform_ODE, np.array([0, 0, 0]), 0)
print(curvedSpring.dBeta)

print(err)
assert(np.isclose(lin.norm(err),0))

######################## Straight Beam Unit Test #######################
# Compare the deflections of a straight, quasi-straight, and           #
# theoretical beam under pure bending moment and ensure they agree     #
########################################################################

print("#################### UNIT TEST 1: Checking Math in straight and quasi-straight cases ####################")
# this test torque is arbitrary
testTorque = 5000

# perform a single pass of solving the deformation ODE under a pure moment
err, res = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([0,0,testTorque]),testTorque)
err, res = quasiStraightSpring.forward_integration(quasiStraightSpring.deform_ODE, np.array([0,0,testTorque]),testTorque)

# plot the results of bending the straight spring to ensure expected result
straightSpring.full_results()
# get the geometry results for the quasi straight beam as well
quasiStraightSpring.full_results(0,1)

# get resultant displacements to compare:
straightResDisp = np.sqrt((straightSpring.deformedNeutralSurface[-1,-1]-straightSpring.undeformedNeutralSurface[-1,-1])**2+
                          (straightSpring.deformedNeutralSurface[0,-1]-straightSpring.undeformedNeutralSurface[0,-1])**2)
print("tip displacement for straight beam:               ",straightResDisp)
# evaluaeory sote beam thlution:
radiusCurvature = (straightSpring.E*straightSpring.IcPts[0]/testTorque)
theoryDisplacement = radiusCurvature-np.sqrt(radiusCurvature**2-straightSpring.fullParamLength**2)
print("tip displacement for straight beam by beam theory:",theoryDisplacement)
# compare with quasi straight beam
quasiStraightResDisp = np.sqrt((quasiStraightSpring.deformedNeutralSurface[-1,-1]-quasiStraightSpring.undeformedNeutralSurface[-1,-1])**2+
                               (quasiStraightSpring.deformedNeutralSurface[0,-1]-quasiStraightSpring.undeformedNeutralSurface[0,-1])**2)
print("tip displacement for quasi-straight beam:         ",quasiStraightResDisp)
print("the above three values should be close to the same")

######################## Curved Beam Unit Test #######################
# Test the deformation functions all together, and use either method #
# to deform a beam under a target torque load.                       #
######################################################################

print("#################### UNIT TEST 2: Deforming Single Spring ####################")
# this torque is arbitrary but in the ballpark for capacity on a gen2 size 6
# actuator
testTorque = 5000
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

    curvedSpring.generate_undeformed_surfaces()
    curvedSpring.calculate_stresses()

    # compare the two guesses
    print("solution planar force vector:", SF)
    # print out results
    print("torque:", testTorque)
    print("angular deformation:", curvedSpring.dBeta/deg2rad)
    print("stiffness (lbf/deg)", testTorque/(curvedSpring.dBeta/deg2rad))
    print("max stress:", curvedSpring.maxStress)
    print("des stress:", curvedSpring.designStress)
    # if it didnt diverge consider the deform ation succesful
    # deformBool = not(divergeFlag)
    deformBool = True
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
np.savetxt("surfaces\A_surface.txt", A, delimiter=",", fmt='%f')
np.savetxt("surfaces\B_surface.txt", B, delimiter=",", fmt='%f')

A = np.hstack((straightSpring.undeformedASurface,np.atleast_2d(np.zeros(len(straightSpring.undeformedASurface))).T))
B = np.hstack((straightSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(straightSpring.undeformedBSurface))).T))
np.savetxt("surfaces\A_surface_straight.txt", A, delimiter=",", fmt='%f')
np.savetxt("surfaces\B_surface_straight.txt", B, delimiter=",", fmt='%f')

##### This section is not currently used but visualizes the "smart initial guess" linear extrapolation process

# make_SW_part(funnySpring)
# print("called SW function")
## Evaluate the curve of angles and magnitudes ###

# print("### TESTING SMART GUESS RULES ###")
# for i in np.linspace(5,testTorque,5):
#     garbage, SF, trash, extra = curvedSpring.wrapped_torque_deform(i,curvedSpring.deform_ODE)
#     print(SF)
# angles = np.empty(len(curvedSpring.wrapped_torque_deform.all_output))
# magnitudes = np.empty(len(curvedSpring.wrapped_torque_deform.all_output))
# for i in range(len(curvedSpring.wrapped_torque_deform.all_output)):
#     angles[i] = np.arctan2(curvedSpring.wrapped_torque_deform.all_output[i][1],curvedSpring.wrapped_torque_deform.all_output[i][0])
#     magnitudes[i] = lin.norm(curvedSpring.wrapped_torque_deform.all_output[i][0:2])
# plt.figure(999)
# plt.scatter(curvedSpring.wrapped_torque_deform.all_moments, angles)
# plt.scatter(curvedSpring.wrapped_torque_deform.all_moments, magnitudes)
# plt.axhline(curvedSpring.betaAngles[-1])
# plt.axhline(np.average(angles))

# print(np.hstack((np.atleast_2d(np.linspace(5,5000,101)).T,np.atleast_2d(magnitudes).T,np.atleast_2d(angles).T)))
plt.show()
