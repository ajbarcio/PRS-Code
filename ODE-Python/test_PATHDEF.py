import pytest
import inspect
from PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control

import numpy as np
from numpy import linalg as lin

from utils import deg2rad

"""
This test file tests all of the currently active path definitions you can choose from. To make your own test, definitely include the last line of each test
to make sure you have run and can call all of the common (i.e., abstractmethod) methods and access (therefore) all of the common attributes
"""

# @pytest.mark.skip(reason="This test is temporarily skipped")
# @pytest.mark.order(1)
def test_definingPoints():
    testPath = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3))
    assert testPath.initialRadius == lin.norm((1,1))
    assert testPath.finalRadius   == lin.norm((-1,3))
    assert testPath.arcLen is not None
    assert np.all([np.isclose(a,b) for a, b in zip(testPath.endPoint, (testPath.get_xy_n(testPath.arcLen,'x'),testPath.get_xy_n(testPath.arcLen,'y')))])
    print(testPath.measure_length(), end=" ")
    # Run all abstract methods in the path under test
    [method() for method in testPath.__class__.__abstractmethods__ if callable(getattr(testPath, method))]

# @pytest.mark.skip(reason="This test is temporarily skipped")
# @pytest.mark.order(2)
def test_definingRadii():
    testPath = LinearRnSpiral(2, np.pi/2, initialRadius=1, finalRadius=3)
    assert testPath.startPoint==(testPath.initialRadius,0)
    assert testPath.endPoint==(0,testPath.finalRadius)
    assert testPath.arcLen is not None
    assert np.all([np.isclose(a,b) for a, b in zip(testPath.endPoint, (testPath.get_xy_n(testPath.arcLen,'x'),testPath.get_xy_n(testPath.arcLen,'y')))])
    assert True
    print(testPath.measure_length(), end=" ")
    # Run all abstract methods in the path under test
    [method() for method in testPath.__class__.__abstractmethods__ if callable(getattr(testPath, method))]

# @pytest.mark.skip(reason="This test is temporarily skipped")
# @pytest.mark.order(3)
def test_definingAll():
    testPath = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3), initialRadius=1, finalRadius=3)
    assert testPath.startPoint == (1,1)
    assert testPath.endPoint == (-1,3)
    assert testPath.initialRadius == 1
    assert testPath.finalRadius == 3
    assert testPath.arcLen is not None
    # print("\n")
    assert(np.isclose(lin.norm(np.subtract(testPath.centerOfCurvature,testPath.startPoint)),testPath.initialRadius))
    assert(np.isclose(lin.norm(np.subtract(testPath.centerOfCurvature,testPath.endPoint)),testPath.finalRadius))
    assert np.all([np.isclose(a,b) for a, b in zip(testPath.startPoint, (testPath.get_xy_n(0,'x'),testPath.get_xy_n(0,'y')))])
    assert np.all([np.isclose(a,b) for a, b in zip(testPath.endPoint, (testPath.get_xy_n(testPath.arcLen,'x'),testPath.get_xy_n(testPath.arcLen,'y')))])
    print(testPath.measure_length(), end=" ")
    assert testPath.momentArm == tuple(np.subtract(testPath.endPoint, testPath.startPoint))
    # Run all abstract methods in the path under test
    [method() for method in testPath.__class__.__abstractmethods__ if callable(getattr(testPath, method))]

# @pytest.mark.skip(reason="This test is temporarily skipped")
# @pytest.mark.order(4)
def test_defaultSpring():
    testPath = RadiallyEndedPolynomial(n=2, arcLen=5.911)
    assert(np.isclose(testPath.arcLen,testPath.measure_length()))
    assert np.all([np.isclose(a,b) for a, b in zip(testPath.startPoint, (testPath.get_xy_n(0,'x'),testPath.get_xy_n(0,'y')))])
    assert np.all([np.isclose(a,b) for a, b in zip(testPath.endPoint, (testPath.get_xy_n(testPath.arcLen,'x'),testPath.get_xy_n(testPath.arcLen,'y')))])
    assert not hasattr(testPath.get_rn(testPath.arcLen), "__len__")
    assert hasattr(testPath.get_rn([0,testPath.arcLen]), "__len__")
    print(testPath.measure_length(), end=" ")
    # Run all abstract methods in the path under test
    [method() for method in testPath.__class__.__abstractmethods__ if callable(getattr(testPath, method))]