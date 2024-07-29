import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring_old import Spring
from materials import Maraging300Steel

name = "manual_spring"
path = "springs\\"+name
parameters = np.genfromtxt(path, delimiter=',')
# print(parameters)
n = parameters[0]
fullParamLength = parameters[1]
thick = parameters[2]
radii = parameters[3:7]
betaAngles = parameters[7:11]
IcPts = parameters[11:15]
IcParamLens = parameters[15:17]
XYParamLens = parameters[17:19]
## ^^ this is dumb but who cares
## generate spring starting guess
tunedSpring = Spring(Maraging300Steel(), n=int(n), radii=radii,
                                         betaAngles=betaAngles, IcPts=IcPts,
                                         IcParamLens=IcParamLens,
                                         XYParamLens=XYParamLens, name="tunedSpring")
# hold inner and outer radii constant, hold XY lengths constant, and beta0
discludeVector = np.ones(len(tunedSpring.parameterVector))
discludeVector[0] = 0
discludeVector[3] = 0
discludeVector[7] = 0
discludeVector[-2] = 0
discludeVector[-1] = 0
print(discludeVector)
tunedSpring.tune_stiffness(910.92, 4554.6, discludeVector)
# tunedSpring =