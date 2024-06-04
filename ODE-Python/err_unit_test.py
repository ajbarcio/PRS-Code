import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring import Spring

deg2rad = np.pi/180

################################################################################
# This script tests the new method of finding angles between two vectors and
# its effect on evaluating the error in the boundary conditions
################################################################################

# Fit the Gen 2 actuator form factor
R0 = 2.2/2
R3 = 5.9/2

# All these radii, angles, and Ic come from a previous version of this code and
# can be considered as arbitrary starting points

R1 = (R0+R3)/2+.26+.125+.125
R2 = (R0+R3)/2-.25+.125

fullAngle = 180

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

XYParamLens = np.array([0.333*0.8,0.667*1.2])

# Generate the spring:

curvedSpring = Spring(n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens)

#### Unit Test 1: Ensure deflection angle is 0 for zero force ####
####              Ensure deflection angle is right sign       ####
####                     (use a small torque)                 ####

for i in range(360):
    fullAngle = i

    beta1 = fullAngle/3*deg2rad*.5*1.2

    beta2 = 2*fullAngle/3*deg2rad*0.9*.9
    beta0 = fullAngle*deg2rad
    inputBetaAngles = np.array([0,beta1,beta2,beta0])
    curvedSpring = Spring(n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens)
    curvedSpring.forward_integration(curvedSpring.deform_ODE,
                                 np.array([0,0,0]),
                                 0)
    if not curvedSpring.dBeta==0:
        print(curvedSpring.dBeta)
    assert(curvedSpring.dBeta==0)
    curvedSpring.forward_integration(curvedSpring.deform_ODE,
                                 np.array([0,0,100]),
                                 100)
    assert(np.sign(curvedSpring.dBeta)>0)
    curvedSpring.forward_integration(curvedSpring.deform_ODE,
                                 np.array([0,0,100]),
                                 -100)
    assert(np.sign(curvedSpring.dBeta)<0)
    print(i, end="\r")

print("     ",end="\r")
#### Unit Test 2: Ensure small torques converge for all angles ####

for i in range(360-120):
    fullAngle = (i+120)

    beta1 = fullAngle/3*deg2rad*.5*1.2
    beta2 = 2*fullAngle/3*deg2rad*0.9*.9
    beta0 = fullAngle*deg2rad
    inputBetaAngles = np.array([0,beta1,beta2,beta0])
    # print(inputBetaAngles)
    curvedSpring = Spring(n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens)
    curvedSpring.deform_by_torque(100,curvedSpring.deform_ODE)
    print("                                                             ", end="\r")
    print(curvedSpring.dBeta, i, end="\r")

for i in range(360-120):
    fullAngle = (i+120)

    beta1 = fullAngle/3*deg2rad*.5*1.2
    beta2 = 2*fullAngle/3*deg2rad*0.9*.9
    beta0 = fullAngle*deg2rad
    inputBetaAngles = np.array([0,beta1,beta2,beta0])
    # print(inputBetaAngles)
    curvedSpring = Spring(n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens)
    curvedSpring.deform_by_torque(-100,curvedSpring.deform_ODE)
    print("                                                             ", end="\r")
    print(curvedSpring.dBeta, i, end="\r")

print("                                                             ", end="\r")
print("success \a")