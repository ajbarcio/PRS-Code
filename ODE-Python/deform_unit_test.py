from stiffness_library import *

R0 = 2.2/2
R3 = 5.9/2

R1 = (R0+R3)/2+.26
R2 = (R0+R3)/2-.25

funnyAngle = 150

betaB = funnyAngle/3*deg2rad*.5
betaC = 2*funnyAngle/3*deg2rad
betaD = funnyAngle*deg2rad
beta0 = betaD

Ics = np.array([0.006, 0.0000025, 0.001])

curvedSpring = Spring(n = 1, radii=np.array([R0,R1,R2,R3]),betaAngles=np.array([0,betaB,betaC,beta0]),IcPts=Ics,IcArcLens=np.array([0.5]))

straightSpring = Spring(n = 1, fullArcLength = 6, radii = np.array([1,3,5,7]),
                        betaAngles=np.array([0,1,2,3])*deg2rad,
                        IcPts = np.array([0.03125, 0.03125, 0.03125]), resolution = 200)

## True Unit Test: Make sure that for the degenerate case (straight beam and arbit. beam with 0 torque)
#  the expected values are returned

straightSpring.spring_geometry(plotBool=0,deformBool=0)
err, res = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([0,0,0]), 0)
# print(err)
assert(np.isclose(lin.norm(err),0))
curvedSpring.spring_geometry(plotBool=0,deformBool=0)
err, res = curvedSpring.forward_integration(curvedSpring.deform_ODE, np.array([0,0,0]), 0)
# print(err)
assert(np.isclose(lin.norm(err),0))

## Straight Beam Unit Test: see if the straight beam behavior matches theory

# print(straightSpring.dxids)
testTorque = 1000
err, res = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([0,0,testTorque]),testTorque)
straightSpring.spring_geometry(plotBool=1,deformBool=1)
print(straightSpring.deformedNeutralSurface[-1,-1]-straightSpring.undeformedNeutralSurface[-1,-1])
r = (straightSpring.E*0.03125/testTorque)
print(r-np.sqrt(r**2-straightSpring.fullArcLength**2))
# assert(np.isclose(straightSpring.deformedNeutralSurface[-1,-1],r-np.sqrt(r**2-straightSpring.fullArcLength**2)))

## Not really a unit test, but still a test
testTorque = 2500
method = "smartGuess"
if method=="smartGuess":
    SFGuess = curvedSpring.smart_initial_load_guess(testTorque,curvedSpring.deform_ODE)
    print(SFGuess)

    res, SF, divergeFlag, i = curvedSpring.wrapped_torque_deform(testTorque,curvedSpring.deform_ODE,SF=SFGuess)
    print(SF)

    deformBool = not(divergeFlag)
    print(curvedSpring.dBeta/deg2rad)
    curvedSpring.spring_geometry(deformBool=deformBool, plotBool=True)
elif method=="slowRamp":
    res, SF, divergeFlag, i = curvedSpring.deform_by_torque_slowRamp(testTorque,curvedSpring.deform_ODE)
    deformBool = not(divergeFlag)
    curvedSpring.spring_geometry(deformBool=deformBool)

A = np.hstack((curvedSpring.undeformedASurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedASurface))).T))
B = np.hstack((curvedSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedBSurface))).T))
np.savetxt("A_surface.txt", A, delimiter=",", fmt='%f')
np.savetxt("B_surface.txt", B, delimiter=",", fmt='%f')

A = np.hstack((straightSpring.undeformedASurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedASurface))).T))
B = np.hstack((straightSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(curvedSpring.undeformedBSurface))).T))
np.savetxt("A_surface_straight.txt", A, delimiter=",", fmt='%f')
np.savetxt("B_surface_straight.txt", B, delimiter=",", fmt='%f')



# make_SW_part(funnySpring)
# print("called SW function")

### Evaluate the curve of angles and magnitudes ###
# for i in np.linspace(5,testTorque,5):
#     fuck, SF, you, dad = curvedSpring.wrapped_torque_deform(i,curvedSpring.deform_ODE)
#     print(SF)
# angles = np.empty(len(curvedSpring.wrapped_torque_deform.all_output))
# magnitudes = np.empty(len(curvedSpring.wrapped_torque_deform.all_output))
# for i in range(len(curvedSpring.wrapped_torque_deform.all_output)):
#     angles[i] = np.arctan2(curvedSpring.wrapped_torque_deform.all_output[i][1],curvedSpring.wrapped_torque_deform.all_output[i][0])
#     magnitudes[i] = lin.norm(curvedSpring.wrapped_torque_deform.all_output[i][0:2])
# plt.figure(999)
# plt.scatter(curvedSpring.wrapped_torque_deform.all_moments, angles)
# plt.scatter(curvedSpring.wrapped_torque_deform.all_moments, magnitudes)

# print(np.hstack((np.atleast_2d(np.linspace(5,5000,101)).T,np.atleast_2d(magnitudes).T,np.atleast_2d(angles).T)))
plt.show()
