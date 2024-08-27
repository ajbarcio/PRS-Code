import numpy as np
from numpy import linalg as lin
import sympy.matrices as Lin
from matplotlib import pyplot as plt
import sympy as sp

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

from sympy.solvers.ode.systems import dsolve_system

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
testCrsc1 = CRSCDEF.Piecewise_Ic_Control(pathDef=testPath,
                       IcPts = np.array([.0003, 0.000026, .0003]),
                       IcParamLens = np.array([0.5]))

# fuck this
testPath.get_crscRef(testCrsc1)

# Initialize a spring and generate its shape (rootfinding method)
testSpring = Spring(testCrsc1, materials.Maraging300Steel, 
                    resolution=1000, name="20270730_spring")

t = testSpring.t

print(testSpring.path.XCoeffs)
print(testSpring.path.YCoeffs)
print(testSpring.crsc.IcCoeffs)
print(testSpring.crsc.domains)

XCoeffs = testSpring.path.XCoeffs.flatten()
YCoeffs = testSpring.path.YCoeffs.flatten()
IcCoeffs = testSpring.crsc.IcCoeffs.flatten()

s = sp.symbols('s')
la = sp.Function('la')(s)
lb = sp.Function('lb')(s)
yPoly = sp.Poly(YCoeffs, s).as_expr()
xPoly = sp.Poly(XCoeffs, s).as_expr()
Ic    = sp.Poly(IcCoeffs, s).as_expr()

dyds   = sp.diff(yPoly, s)
dxds   = sp.diff(xPoly, s)
d2yds2 = sp.diff(yPoly, s, 2)
d2xds2 = sp.diff(xPoly, s, 2)

dIcds  = sp.diff(Ic   , s)

rn = ((1+dyds**2/dxds**2)/(d2yds2/dxds-d2xds2*dyds/dxds**2))
drnds = sp.diff(rn, s)
print("-------------")
p_dot = Lin.Matrix([[dIcds],[drnds]])

# f geofunc
f1 = Lin.Matrix([[1/(rn*t),-Ic/(rn**2*t)]])
g1 = Lin.Matrix([[0,1/(rn-la)-1/(rn+lb)-(la+lb)/rn**2]])
p_hat = f1.col_join(g1)

q_hat = Lin.Matrix([[-la, lb],[1/(rn-la)-1/rn, 1/(rn+la)-1/rn]])
# print((q_hat))

LHS = p_hat*p_dot
print(sp.shape(LHS))
RHS = q_hat.LUsolve(LHS)
print(sp.shape(RHS))

la_eqn = sp.diff(la, s)-RHS[0]
lb_eqn = sp.diff(lb, s)-RHS[1]

dsolve_system([la_eqn, lb_eqn])