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
IcLens=np.array([0.4, 0.667])

XYParamLens = np.array([0.333,0.667])

# Generate the spring:

curvedSpring = Spring(n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens)

testTorque = 4554.6

# study sensitivity within 10% of the first radius
print("testing r1, original value:", curvedSpring.parameterVector[1])
results = curvedSpring.sensitivity_study(1,.1,10,4554.6)
np.savetxt("sensitivity_results\sens_to_r1.txt", results, delimiter=",", fmt='%f')

print("testing r2, original value:", curvedSpring.parameterVector[2])
results = curvedSpring.sensitivity_study(2,.1,10,4554.6)
np.savetxt("sensitivity_results\sens_to_r2.txt", results, delimiter=",", fmt='%f')

print("testing b1, original value:", curvedSpring.parameterVector[5])
results = curvedSpring.sensitivity_study(5,.1,10,4554.6)
np.savetxt("sensitivity_results\sens_to_b1.txt", results, delimiter=",", fmt='%f')

print("testing b2, original value:", curvedSpring.parameterVector[6])
results = curvedSpring.sensitivity_study(6,.1,10,4554.6)
np.savetxt("sensitivity_results\sens_to_b2.txt", results, delimiter=",", fmt='%f')

print("testing b3, original value:", curvedSpring.parameterVector[7])
results = curvedSpring.sensitivity_study(1,.1,10,4554.6)
np.savetxt("sensitivity_results\sens_to_b3.txt", results, delimiter=",", fmt='%f')