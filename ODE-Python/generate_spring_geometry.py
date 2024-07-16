import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from StatProfiler import SSProfile
import scipy

from spring import Spring
import spring



# Parameters: (p = number of control points)
#             (q = number of Ic profile points (you don't need these))
#
#         radii: vector of p controlled radii
#    betaAngles: vector of p controlled angles from the horizontal (RADIANS)
#   XYParamLens: vector of p-2 fractions of overall parameter length at which to
#                control neutral surface
#         IcPts: vector of q controlled Ic
#   IcParamLens: vector of q-2 fractions of overall parameter length at which to
#              control Ic
#
# Other Parameters:
#
#                   n = number of "legs"
#     fullParamLength = fake, don't use this (leave at default of 6)
#   outPlaneThickness = thickness of spring (.375 is current Gen2 thickness)
#          resolution = how many points to generate along each curve


# You can use other methods for generating/reading in the parameters, but:
Rinner = 1
Router = 2.5 # these are the form factor for the gen2 actuator
# This is just a placeholder value
Rintermediate = (Rinner+Router)/2

radii = [Rinner, Rintermediate, Rintermediate, Router]

beta0 = 160*spring.deg2rad

betaAngles = [0,beta0*1/3,beta0*2/3,beta0]

# Leave these unchanged if you want to only extract the neutral radius path
IcPts = [0.0025,0.0025]
IcParamLens = []

abaqusSpring = spring()

