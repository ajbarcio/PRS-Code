import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from spring_old import Spring
from materials import Maraging300Steel

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

straightSpring = Spring(Maraging300Steel(), n=6, radii=inputRadii,
                                        betaAngles=np.array([0,0])*deg2rad,
                                        outPlaneThickness=outThick,
                                        IcPts=np.array([IcStraight,IcStraight,IcStraight]),
                                        IcParamLens = np.array([0.5]),
                                        XYParamLens=np.array([1/3.0,2/3.0]),
                                        name="testStraightSpring")

###### TEST 0: Just Forward Integration
print("TEST 0")
testTorque = 4554.6

err0, res0 = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([0,0,testTorque]),testTorque)
straightSpring.full_results()
plt.show()

###### TEST 0.5: Forward Integration with Tension:
print("TEST 0.5")
err1, res1 = straightSpring.forward_integration(straightSpring.deform_ODE, np.array([100,100,testTorque]),testTorque)
straightSpring.full_results()
print(err0)
print(err1-err0)
print(lin.norm(res1-res0))
plt.show()

if False:
####### TEST 1: Just Start from 0
    print("TEST 1")
    res, SF, divergeFlag, i = straightSpring.deform_by_torque(testTorque,straightSpring.deform_ODE)
    straightSpring.full_results()
    plt.show()
