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
res, SF, divergeFlag = funnySpring.deform_by_torque(5000,funnySpring.deform_ODE)
print(SF)
