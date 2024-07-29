import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring_old import Spring
from materials import Maraging300Steel
from StatProfiler import SSProfile

R0 = 2.0/2
R3 = 5.0/2

deg2rad = np.pi/180

# R1 = R3/3
# R2 = 2*R3/3
inputRadii      = np.array([R0,R3])
straightThickness = 0.25
outThick = 0.375
IcStraight = outThick*straightThickness**3/12
# print(IcStraight)

straightSpring = Spring(Maraging300Steel(), n=2, radii=inputRadii,
                                        betaAngles=np.array([0,0])*deg2rad,
                                        outPlaneThickness=outThick,
                                        IcPts=np.array([IcStraight,IcStraight]),
                                        IcParamLens = np.array([]),
                                        XYParamLens=np.array([]),
                                        name="testStraightSpring")

testTorque=4556.4

# Just forward integration
straightSpring.generate_undeformed_surfaces()

err, res0 = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([0,0,0]), 0)
err, res1 = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([100,0,0]), 0)

print(res0)
print(res0[:,-1])
print(res1[:,-1])

