import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring import Spring

deg2rad = np.pi/180

# Fit the Gen 2 actuator form factor
R0 = 2.2/2
R3 = 5.9/2

# All these radii, angles, and Ic come from a previous version of this code and
# can be considered as arbitrary starting points

R1  = (R0+R3)/2
R15 = (R0+R3)/2+.25
R2  = (R0+R3)/2

fullAngle = 150

beta1  = fullAngle/4*deg2rad*.5
beta15 = fullAngle/2*deg2rad*.8
beta2  = 3*fullAngle/4*deg2rad*1.1
beta0  = fullAngle*deg2rad

Ics    = np.array([0.015, 0.0001, 0.008])

curvedSpring = Spring(n = 2, 
                      radii=np.array([R0,R1,R15,R2,R3]),
                      betaAngles=np.array([0,beta1,beta15,beta2,beta0]),
                      IcPts=Ics)
curvedSpring.vis_results(deformBool=0)
