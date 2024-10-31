import pytest
import inspect
from PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from CRSCDEF import IcDef, Constant_Ic, Piecewise_Ic_Control

import numpy as np
from numpy import linalg as lin

from spring import deg2rad

"""
This test file tests all of the currently active path definitions you can choose from. To make your own test, definitely include the last line of each test
to make sure you have run and can call all of the common (i.e., abstractmethod) methods and access (therefore) all of the common attributes
"""

# @pytest.mark.skip(reason="This test is temporarily skipped")
def test_Constant_Ic():
    testPath1 = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3))
    testCrsc  = Constant_Ic(testPath1, 0.375, 0.125)
    [method() for method in testCrsc.__class__.__abstractmethods__ if callable(getattr(testCrsc, method))]

def test_Constant_Ic():
    testPath1 = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3))
    testCrsc  = Constant_Ic(testPath1, 0.375, 0.125)
    [method() for method in testCrsc.__class__.__abstractmethods__ if callable(getattr(testCrsc, method))]
    
# @pytest.mark.skip(reason="This test is temporarily skipped")
def test_definingRadii():
    testPath2 = LinearRnSpiral(2, np.pi/2, initialRadius=1, finalRadius=3)
    assert testPath2.startPoint==(testPath2.initialRadius,0)
    assert testPath2.endPoint==(0,testPath2.finalRadius)
    assert testPath2.arcLen is not None
    assert np.all([np.isclose(a,b) for a, b in zip(testPath2.endPoint, (testPath2.get_xy_n(testPath2.arcLen,'x'),testPath2.get_xy_n(testPath2.arcLen,'y')))])
    assert True
    print(testPath2.measure_length())
    # Run all abstract methods in the path under test
    [method() for method in testPath2.__class__.__abstractmethods__ if callable(getattr(testPath2, method))]

# @pytest.mark.skip(reason="This test is temporarily skipped")
def test_definingAll():
    testPath3 = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3), initialRadius=1, finalRadius=3)
    assert testPath3.startPoint == (1,1)
    assert testPath3.endPoint == (-1,3)
    assert testPath3.initialRadius == 1
    assert testPath3.finalRadius == 3
    assert testPath3.arcLen is not None
    # print("\n")
    assert(np.isclose(lin.norm(np.subtract(testPath3.centerOfCurvature,testPath3.startPoint)),testPath3.initialRadius))
    assert(np.isclose(lin.norm(np.subtract(testPath3.centerOfCurvature,testPath3.endPoint)),testPath3.finalRadius))
    assert np.all([np.isclose(a,b) for a, b in zip(testPath3.startPoint, (testPath3.get_xy_n(0,'x'),testPath3.get_xy_n(0,'y')))])
    assert np.all([np.isclose(a,b) for a, b in zip(testPath3.endPoint, (testPath3.get_xy_n(testPath3.arcLen,'x'),testPath3.get_xy_n(testPath3.arcLen,'y')))])
    print(testPath3.measure_length())
    assert testPath3.momentArm == tuple(np.subtract(testPath3.endPoint, testPath3.startPoint))
    # Run all abstract methods in the path under test
    [method() for method in testPath3.__class__.__abstractmethods__ if callable(getattr(testPath3, method))]

def test_defaultSpring():
    testPath4 = RadiallyEndedPolynomial(n=2, arcLen=5.911)
    assert(np.isclose(testPath4.arcLen,testPath4.measure_length()))
    assert np.all([np.isclose(a,b) for a, b in zip(testPath4.startPoint, (testPath4.get_xy_n(0,'x'),testPath4.get_xy_n(0,'y')))])
    assert np.all([np.isclose(a,b) for a, b in zip(testPath4.endPoint, (testPath4.get_xy_n(testPath4.arcLen,'x'),testPath4.get_xy_n(testPath4.arcLen,'y')))])
    assert not hasattr(testPath4.get_rn(testPath4.arcLen), "__len__")
    assert hasattr(testPath4.get_rn([0,testPath4.arcLen]), "__len__")
    print(testPath4.measure_length())
    # Run all abstract methods in the path under test
    [method() for method in testPath4.__class__.__abstractmethods__ if callable(getattr(testPath4, method))]