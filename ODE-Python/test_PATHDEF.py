import pytest
from PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial

import numpy as np
from numpy import linalg as lin

from spring import deg2rad

@pytest.mark.skip(reason="This test is temporarily skipped")
def test_definingPoints():
    testPath1 = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3))
    assert testPath1.initialRadius == lin.norm((1,1))
    assert testPath1.finalRadius   == lin.norm((-1,3))
    assert testPath1.arcLen is not None
    assert np.all([np.isclose(a,b) for a, b in zip(testPath1.endPoint, (testPath1.get_xy_n(testPath1.arcLen,'x'),testPath1.get_xy_n(testPath1.arcLen,'y')))])
    print("\n", testPath1.measure_length())

@pytest.mark.skip(reason="This test is temporarily skipped")
def test_definingRadii():
    testPath2 = LinearRnSpiral(2, np.pi/2, initialRadius=1, finalRadius=3)
    assert testPath2.startPoint==(testPath2.initialRadius,0)
    assert testPath2.endPoint==(0,testPath2.finalRadius)
    assert testPath2.arcLen is not None
    assert np.all([np.isclose(a,b) for a, b in zip(testPath2.endPoint, (testPath2.get_xy_n(testPath2.arcLen,'x'),testPath2.get_xy_n(testPath2.arcLen,'y')))])
    assert True
    print(testPath2.measure_length())

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

def test_defaultSpring():
    testPath4 = RadiallyEndedPolynomial(n=2, arcLen=6)