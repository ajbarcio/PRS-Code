import pytest
import inspect

import numpy as np
from numpy import linalg as lin
from matplotlib import pyplot as plt

from PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control
from materials import TestMaterial
from spring import Spring

from utils import deg2rad

from StatProfiler import SSProfile

def test_classFunctionality():

    """ This is one huge test to evaluate the overall functionality of the class. This includes:
            proper inheritance form path, crsc, and material objects
            geometry plotting/calcualtion
            deforming a spring of known geometry
                using each of the different methods to attempt to avoid divergence
            plotting a deformed spring
    """

    print("\n")
    # testPath = LinearRnSpiral(2, None, (1,0),(-1.5,2),5,3)
    # testPath = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3), initialRadius=np.sqrt(2), finalRadius=5)
    testPath = RadiallyEndedPolynomial(2,6)
    testCrsc = Constant_Ic(testPath, 0.375, Ic0=0.00125)
    testSprg = Spring(testPath, testCrsc, TestMaterial)
    # Test that the plotting works; you, unfortunately, must be smart about this: if it looks right, success
    testSprg.plot_spring()
    
    
    # Test the trivial deformation case: deform by 0 torque
    testSprg.deform_by_torque(0,testSprg.deform_ODE)
    assert np.isclose(testSprg.dBeta,0)
    
    # Test a minor deformation case: deform by low (5-15NM) torque, ensure no divergence
    res, solnSF, divergeFlag, i  = testSprg.deform_by_torque(5,testSprg.deform_ODE,SF=np.array([0,0,5]))
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSF, 5)
    assert np.isclose(lin.norm(err),0)
    assert divergeFlag == 0
    
    # Test if you can simply deform the spring straight-up:
    SSProfile("nonlinear").tic()
    res, solnSFOrig, divergeFlag, i = testSprg.deform_by_torque(testSprg.torqueCapacity,testSprg.deform_ODE)
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSFOrig, testSprg.torqueCapacity)
    # print(err)
    assert np.isclose(lin.norm(err),0)
    SSProfile("nonlinear").toc()
    print("Starting from 0 takes:", SSProfile("nonlinear").agg, "s")

    SSProfile("predict forces").tic()
    # Test if you can deform the spring by assuming linearity:
    res, solnSF, divergeFlag, i = testSprg.deform_by_torque_predict_forces(testSprg.torqueCapacity,testSprg.deform_ODE)
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSF, testSprg.torqueCapacity)
    if not np.isclose(lin.norm(err),0):
        print("PREDICT FORCES DIDN'T WORK")
    else:
        print("PREDICT FORCES WORKED")
    print(solnSF)
    print(solnSFOrig)
    SSProfile("predict forces").toc()
    print(SSProfile("predict forces").agg)
    print("Trying to predict individual forces takes:", SSProfile("predict forces").agg, "s")

    SSProfile("predict angle").tic()
    res, solnSF, divergeFlag, i = testSprg.deform_by_torque_predict_angle(testSprg.torqueCapacity,testSprg.deform_ODE)
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSF, testSprg.torqueCapacity)
    if not np.isclose(lin.norm(err),0):
        print("PREDICT ANGLE DIDN'T WORK")
    else:
        print("PREDICT ANGLE WORKED")
    print(solnSF)
    print(solnSFOrig)
    SSProfile("predict angle").toc()
    print("Trying to predict force angle takes:", SSProfile("predict angle").agg, "s")

    testSprg.plot_deform(showBool=False)


    # [method() for name, method in testSprg.__class__.__dict__.items() if callable(method) and not name.startswith('__') and not name.startswith('_') for method in [getattr(testSprg, name)]]
    plt.show()
