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

funnySpring = Spring(radii=np.array([R0,R1,R2,R3]),betaAngles=np.array([0,betaB,betaC,beta0]))
# SFGuess = funnySpring.smart_initial_load_guess(50000,funnySpring.deform_ODE)
# print(SFGuess)

# funnySpring.wrapped_torque_deform(50000,funnySpring.deform_ODE,SF=SFGuess)

for i in np.linspace(5,5000,5):
    fuck, SF, you = funnySpring.wrapped_torque_deform(i,funnySpring.deform_ODE)
    print(SF)
angles = np.empty(len(funnySpring.wrapped_torque_deform.all_output))
magnitudes = np.empty(len(funnySpring.wrapped_torque_deform.all_output))
for i in range(len(funnySpring.wrapped_torque_deform.all_output)):
    angles[i] = np.arctan2(funnySpring.wrapped_torque_deform.all_output[i][1],funnySpring.wrapped_torque_deform.all_output[i][0])
    magnitudes[i] = lin.norm(funnySpring.wrapped_torque_deform.all_output[i][0:2])
plt.plot(angles)
plt.plot(magnitudes)

# print(np.hstack((np.atleast_2d(np.linspace(5,5000,101)).T,np.atleast_2d(magnitudes).T,np.atleast_2d(angles).T)))
plt.show()
