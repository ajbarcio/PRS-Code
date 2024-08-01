import numpy as np
from scipy import optimize as opt

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

# Standard procedure: define a path, then define a thickness profile using the
#                     length of that path, then pass both of those objects into
#                     the overall spring object that contains deformation
#                     methods

IR = 1.3
OR = 2.013

testPath = PATHDEF.Minimal_Polynomial_Definition4(n=2, fullParamLength=4,
                                       radii = np.array([IR,(IR+OR)/2*1.15,OR]),
                                       ffradii = np.array([1.14, 2.113]),
                                       alphaAngles = np.array([45,45])*deg2rad,
                                       betaAngles = np.array([0,87.5,175])*deg2rad,
                                       XYFactors = np.array([0.5]))
testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = np.array([.002, .000026, .000026, .0003]),
                       IcParamLens = np.array([.50, .6]))
# fuck this
testPath.get_crscRef(testCrsc)

testSpring = Spring(testCrsc, materials.Maraging300Steel, name="20270730_spring")

# Actually try to deform it
# res, SF, divergeFlag, i = testSpring.deform_by_torque(-4549,
#                                                       testSpring.deform_ODE,
#                                                       SF=np.array([0,0,-4549]))
# testSpring.plot_deform(showBool=False)
res, SF, divergeFlag, i = testSpring.deform_by_torque(4549,
                                                      testSpring.deform_ODE,
                                                      SF=np.array([0,0,4549]))

print(testSpring.solnerr)
print(SF)
testSpring.plot_deform()

testSpring.export_surfaces()

#just try to print it out
# testSpring.plot_spring(showBool=True)

# def Ic_optimizer(IcParams):
#     lowerIcParam = IcParams[0]
#     upperIcParam = IcParams[0]+IcParams[1]
#     testCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
#                        IcPts = np.array([.002, .00004, .00004, .002*.5]),
#                        IcParamLens = np.array([lowerIcParam, upperIcParam]))
#     testPath.get_crscRef(testCrsc)
#     testSpring = Spring(testCrsc, materials.Maraging300Steel, name="test_spring")
#     deflection = stiffness_eval(testSpring,4549)
#     return 1/deflection


# def stiffness_eval(spring, testTorque):
#     spring.deform_by_torque(-testTorque, spring.deform_ODE,SF=np.array([0,0,-testTorque]))
#     dBeta1 = spring.dBeta
#     spring.deform_by_torque(testTorque, spring.deform_ODE,SF=np.array([0,0,testTorque]))
#     dBeta2 = spring.dBeta
#     dBetaRes = np.max([dBeta1,dBeta2])
#     return dBetaRes

# print("about to try to optimize")
# res = opt.minimize(Ic_optimizer, [.5, .1])
# print(res)
