import pytest
import inspect
from PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control

import numpy as np
from numpy import linalg as lin

from spring import deg2rad

# @pytest.mark.skip(reason="This test is temporarily skipped")
# @pytest.mark.order(5)
def test_Constant_Ic():
    # Make a crsc profile for some simple curved path
    testPath = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3))
    testCrsc  = Constant_Ic(testPath, 0.375, Ic0=0.125)

    # Spot test your path tho ensure your math is accurate and behaving
    check = checkMath(testCrsc, testPath)
    assert check
    # Run all the standard functions lmao
    [method() for method in testCrsc.__class__.__abstractmethods__ if callable(getattr(testCrsc, method))]

# @pytest.mark.skip(reason="This test is temporarily skipped")
# @pytest.mark.order(6)
def test_Piecewise_Ic():
    testPath = RadiallyEndedPolynomial(2, 6)
    testCrsc = Piecewise_Ic_Control(testPath, 0.375, IcPts=np.array([.005,.0005,.005]),IcParamLens=np.array([0.5]))

    # Spot test your path tho ensure your math is accurate and behaving
    check = checkMath(testCrsc, testPath)
    assert check
    # Run all the standard functions lmao

    [method() for method in testCrsc.__class__.__abstractmethods__ if callable(getattr(testCrsc, method))]
    
def checkMath(crsc, path):
    randomCoord = np.random.rand()*crsc.arcLen
    print("Checking compliance at xi = ", randomCoord, end=" ")
    lalb = crsc.get_lalb(randomCoord)
    la = lalb[0]
    lb = lalb[1]
    Ic   = crsc.get_Ic(randomCoord)
    e    = crsc.get_eccentricity(randomCoord)
    h    = crsc.get_Thk(randomCoord)
    rn   = path.get_rn(randomCoord)
    # These two equalities should hold for valid geometry
    check1 =  np.isclose((Ic),(crsc.t*rn/2*(lb**2-la**2)))
    check2 =  np.isclose((rn),((la+lb)/(np.log((rn+lb)/(rn-la)))))
    return check1 and check2