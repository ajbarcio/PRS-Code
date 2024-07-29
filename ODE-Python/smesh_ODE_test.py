import numpy as np
from numpy import linalg as lin
import matplotlib.pyplot as plt
from spring_old import Spring, Deform_Wrapper
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

XYParamLens = np.array([0.333*0.8,0.667*1.2])

# Generate the spring:

testSpring = Spring(Maraging300Steel(), n = 2, radii=inputRadii,
                             betaAngles=inputBetaAngles,
                             IcPts=Ics,
                             IcParamLens=IcLens, name="testSpring")

testSpring.full_results(deformBool=0)
print(testSpring.fullParamLength)
print(testSpring.measure_length())
print(lin.norm(testSpring.dxids-np.ones(len(testSpring.dxids))))

plt.show()
