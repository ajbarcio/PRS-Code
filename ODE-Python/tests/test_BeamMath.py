import pytest
import inspect

import numpy as np
from numpy import linalg as lin
from matplotlib import pyplot as plt

from modules.PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from modules.CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control
from modules.materials import TestMaterial
from modules.spring import Spring, determineFastestSolver

from modules.utils import deg2rad

from modules.StatProfiler import SSProfile

# @pytest.mark.skip()
def test_tensive_straight_beam():
    testPath = RadiallyEndedPolynomial(1, 10, radii=np.array([0, 10]), ffradii=np.array([0, 10]), alphaAngles=np.array([0,0]), betaAngles=np.array([0,0]), XYFactors=np.array([]))
    # testCrsc = Constant_Ic(testPath, 0.375, Ic0=0.00125)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([1/120,1/120]), t=0.1)
    testSprg = Spring(testPath, testCrsc, TestMaterial)


    print("In this plot you should see a deformation of 1 inch")
    err, res = testSprg.forward_integration(testSprg.deform_withTension_ODE, np.array([1000,0,0]), 0)
    
    assert np.isclose(res[1,-1], 11)
    
    testSprg.plot_deform(testSprg)

@pytest.mark.skip()
def test_tensive_curved_beam():
    testPath = RadiallyEndedPolynomial(1, 6, radii = np.array([1, 3]), ffradii = np.array([1, 3]), alphaAngles=np.array([90,90])*deg2rad, betaAngles=np.array([0,180])*deg2rad, XYFactors=np.array([]))
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([1/120,1/120]), t=0.5)
    testSprg = Spring(testPath, testCrsc, TestMaterial)

    testSprg.forward_integration(testSprg.deform_withTension_ODE, np.array([0,10,0]), 0)
    testSprg.plot_deform(testSprg)


# @pytest.mark.skip()
def test_bending_straight_beam():
    testPath = RadiallyEndedPolynomial(1, 10, radii=np.array([0, 10]), ffradii=np.array([0, 10]), alphaAngles=np.array([0,0]), betaAngles=np.array([0,0]), XYFactors=np.array([]))
    # testCrsc = Constant_Ic(testPath, 0.375, Ic0=0.00125)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([0.1,0.1]), t=0.1)
    testSprg = Spring(testPath, testCrsc, TestMaterial, resolution=1000)

    L = (testPath.measure_length())
    P = 30
    E = TestMaterial.E
    I = 0.1

    tipDeflection = P*L**3/(3*E*I)
    print(tipDeflection)
    
    err, res = testSprg.forward_integration(testSprg.deform_ODE, np.array([0,30,0]), 0)

    print(res[2,-1])

    assert np.isclose(res[2,-1], tipDeflection, rtol=0.01)

    # print(tipDeflection*2)
    print("In this plot you should see a tip deflection of 1 inch")
    testSprg.plot_deform(testSprg)

def test_bending_uniform_curved_beam():
    R = 1.4564
    M = 67.81002
    Ic = 0.005
    testPath = LinearRnSpiral(1, R*np.pi, startPoint=(R,0), endPoint=(-R,0),
                                       initialRadius=R,  finalRadius=R)
    testCrsc = Constant_Ic(testPath, 0.375, Ic0 = Ic)
    testSprg = Spring(testPath, testCrsc, TestMaterial, resolution = 200)

    quantity = 1/R+M/TestMaterial.E/Ic
    alpha0 = testPath.get_alpha(0)

    predictedX = np.sin(quantity*np.pi*R+alpha0)/quantity+testSprg.x0-np.sin(alpha0)/quantity
    predictedY = -np.cos(quantity*np.pi*R+alpha0)/quantity-testSprg.y0+np.cos(alpha0)/quantity

    predictDX  = predictedX-testSprg.xL
    predictDY  = predictedY-testSprg.yL
    predictDist = lin.norm((predictDX, predictDY))


    print("")
    print("Theory-predicted displacement")
    print(predictDX, predictDY)
    print(predictDist)

    err, res = testSprg.forward_integration(testSprg.deform_ODE, np.array([0,0,M]), M)

    DX = res[1,-1]-testSprg.xL
    DY = res[2,-1]-testSprg.yL

    print(DX, DY)

    print(testSprg.dist)

    testSprg.plot_deform(False)

    assert(np.isclose(testSprg.dist, predictDist))
    assert(np.isclose(DX, predictDX))
    
