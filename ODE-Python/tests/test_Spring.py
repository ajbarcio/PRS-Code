import pytest
import inspect

import json

import numpy as np
from numpy import linalg as lin
from matplotlib import pyplot as plt
import os

from modules.PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from modules.CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control
from modules.materials import TestMaterial
from modules.spring import Spring, determineFastestSolver

from modules.utils import deg2rad

from modules.StatProfiler import SSProfile

import difflib

# @pytest.mark.skip()
def test_saveAndLoad():
    examplePath = RadiallyEndedPolynomial(2,6)
    exampleCrsc = Piecewise_Ic_Control(examplePath, IcPts=np.array([.008,.0008,.008]))
    exampleSprg = Spring(examplePath, exampleCrsc, TestMaterial, name="Example_Spring")

    beforeAll = exampleSprg.export_parameters()

    loadedSprg = Spring.from_param("springs/Example_Spring")
    loadedSprg.name = loadedSprg.name+"1"
    afterAll = loadedSprg.export_parameters()

    with open("springs/Example_Spring_crsc.param") as before:
        beforeCrsc = json.load(before) 
    with open("springs/Example_Spring1_crsc.param") as after:
        afterCrsc = json.load(after)

    with open("springs/Example_Spring_path.param") as before:
        beforePath = json.load(before)
    with open("springs/Example_Spring1_path.param") as after:
        afterPath = json.load(after)

    with open("springs/Example_Spring_sprg.param") as before:
        beforeSprg = json.load(before)
    with open("springs/Example_Spring1_sprg.param") as after:
        afterSprg = json.load(after)

    for key, value in beforePath.items():
        if key == "arcLen":
            pass
        else:
            assert value==afterPath[key]

    assert beforeCrsc==afterCrsc
    for key, value in beforeSprg.items():
        if key == "name":
            pass
        else:
            assert value==afterSprg[key]

    for a in beforeAll + afterAll:
        try:
            os.remove(a)
        except FileNotFoundError:
            print(f"File {a} not found.")
        except PermissionError:
            print(f"Permission denied to remove {a}.")
        except Exception as e:
            print(f"Error occurred: {e}")
    # assert beforeSprg==afterSprg
    # pass

@pytest.mark.skip(reason="don't need to run this every time (main functionality included in test_classFunctionality)")
def test_modeSelector():
    testPath = RadiallyEndedPolynomial(2,6)
    # testCrsc = Constant_Ic(testPath, 0.375, Ic0=0.00125)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([.008,.0008,.008]))
    testSprg = Spring(testPath, testCrsc, TestMaterial)

    fastestSolver = determineFastestSolver(testSprg)
    print(testSprg.deformMode)
    print(fastestSolver)
    assert fastestSolver==str(testSprg.deformMode)

# @pytest.mark.skip(reason="This test takes eons to run (by that I mean ~12 seconds)")
def test_classFunctionality():

    passBool = 0

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
    # testPath = RadiallyEndedPolynomial(2,6)
    # # testCrsc = Constant_Ic(testPath, 0.375, Ic0=0.00125)
    # testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([.008,.0008,.008]))
    # testSprg = Spring(testPath, testCrsc, TestMaterial)
    # Test that the plotting works; you, unfortunately, must be smart about this: if it looks right, success
    testSprg = Spring.from_param("springs/Size_6_Special_Titanium5")
    testSprg.plot_spring()
    # plt.show()
    
    # Test the trivial deformation case: deform by 0 torque
    testSprg.deform_by_torque(0,testSprg.deform_ODE)
    assert np.isclose(testSprg.dBeta,0)
    print("NO TORQUE WORKED")
    
    # Test a minor deformation case: deform by low (5-15NM) torque, ensure no divergence
    res, solnSF, divergeFlag, i  = testSprg.deform_by_torque(5,testSprg.deform_ODE,SF=np.array([0,0,5]))
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSF, 5)
    assert np.isclose(lin.norm(err),0)
    assert divergeFlag == 0
    print("SMALL TORQUE WORKED")
    # print(res[0,-1])
    # print(testSprg.dBeta)
    # testSprg.plot_deform(showBool=False)
    
    testDeformLevel = testSprg.torqueCapacity # Eventually this needs to work for gain = 1

    SSProfile("everything").tic()
    # Test if you can simply deform the spring straight-up:
    SSProfile("nonlinear").tic()
    res, solnSF, divergeFlag, i = testSprg.deform_by_torque(testDeformLevel,testSprg.deform_ODE)
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSF, testDeformLevel)
    # print(err)
    if not np.isclose(lin.norm(err),0):
        print("NONLINEAR DIDN'T WORK")
    else:
        print("NONLINEAR WORKED")
        testSprg.plot_deform(showBool=False)
        passBool = 1
    print(solnSF)
    print(testSprg.solnerr)
    SSProfile("nonlinear").toc()
    print("Starting from 0 takes:", SSProfile("nonlinear").agg, "s", i, "iterations")

    SSProfile("predict forces").tic()
    # Test if you can deform the spring by assuming linearity:
    res, solnSF, divergeFlag, i = testSprg.deform_by_torque_predict_forces(testDeformLevel,testSprg.deform_ODE,breakBool=False)
    # print(divergeFlag, i)
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSF, testDeformLevel)
    if not np.isclose(lin.norm(err),0):
        print("PREDICT FORCES DIDN'T WORK")
    else:
        print("PREDICT FORCES WORKED")
        testSprg.plot_deform(showBool=False)
        passBool = 1
    print(solnSF)
    print(testSprg.solnerr)
    # print(solnSFOrig)
    SSProfile("predict forces").toc()
    # print(SSProfile("predict forces").agg)
    print("Trying to predict individual forces takes:", SSProfile("predict forces").agg, "s", i+2, "iterations")

    SSProfile("predict angle").tic()
    res, solnSF, divergeFlag, i = testSprg.deform_by_torque_predict_angle(testDeformLevel,testSprg.deform_ODE)
    err, res = testSprg.forward_integration(testSprg.deform_ODE, solnSF, testDeformLevel)
    if not np.isclose(lin.norm(err),0):
        print("PREDICT ANGLE DIDN'T WORK")
    else:
        print("PREDICT ANGLE WORKED")
        testSprg.plot_deform(showBool=False)
        passBool = 1
    print(res[0,-1])
    print(testSprg.dBeta)
    print(solnSF)
    # print(solnSFOrig)
    SSProfile("predict angle").toc()
    print("Trying to predict force angle takes:", SSProfile("predict angle").agg, "s", i+2, "iterations")

    SSProfile("everything").toc()
    print("total time:", SSProfile("everything").agg, "s")
    # testSprg.detailed_deform_regression(testSprg.torqueCapacity,testSprg.deform_ODE,50,0.5,)

    # if not np.sign(res[0,-1])==np.sign(testSprg.dBeta):
    #     raise ValueError("WTF") 
    # [method() for name, method in testSprg.__class__.__dict__.items() if callable(method) and not name.startswith('__') and not name.startswith('_') for method in [getattr(testSprg, name)]]
    plt.show()
    assert passBool

# @pytest.mark.skip(reason="Probably don't skip this one unless you want to only shove springs one way")
def test_negativeTroque():
    pass