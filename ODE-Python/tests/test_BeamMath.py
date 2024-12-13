import pytest
import inspect

import numpy as np
from numpy import linalg as lin
from matplotlib import pyplot as plt

from modules.PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from modules.CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control
from modules.materials import TestMaterial
from modules.spring import Spring, determineFastestSolver

import modules.utils as utils
from modules.utils import deg2rad

from modules.StatProfiler import SSProfile

# @pytest.mark.skip()
def test_tensive_straight_beam():

    """ This test ensures that the tensive aspect of the simulator (optional) agrees with theory """

    # Define a straight beam with uniform second moment of area (therefore, thickness)
    testPath = RadiallyEndedPolynomial(1, 10, radii=np.array([0, 10]), ffradii=np.array([0, 10]), alphaAngles=np.array([0,0]), betaAngles=np.array([0,0]), XYFactors=np.array([]))
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([1/120,1/120]), t=0.1)
    # The material 'TestMaterial' is used to simplify the process of attaining a loading which results in unit displacement
    testSprg = Spring(testPath, testCrsc, TestMaterial)

    # print("In this plot you should see a deformation of 1 inch")
    # Apply a purely tensive force to the beam
    err, res = testSprg.forward_integration(testSprg.deform_withTension_ODE, np.array([1000,0,0]), 0)
    
    # Ensure the beam has elongated by one inch
    assert np.isclose(res[1,-1], 11)
    
    # testSprg.plot_deform(testSprg)

# @pytest.mark.skip()
def test_bending_straight_beam():

    """ This test ensures that a straight beam (no tensive component) bends according to theory under the simulator """

    # Define a straight beam with uniform second moment of area (therefore, thickness)
    testPath = RadiallyEndedPolynomial(1, 10, radii=np.array([0, 10]), ffradii=np.array([0, 10]), alphaAngles=np.array([0,0]), betaAngles=np.array([0,0]), XYFactors=np.array([]))
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([0.1,0.1]), t=0.1)
    testSprg = Spring(testPath, testCrsc, TestMaterial, resolution=1000)

    # Calculate tip deflection indicated by beam theory
    L = (testPath.measure_length())
    P = 30
    E = TestMaterial.E
    I = 0.1

    tipDeflection = P*L**3/(3*E*I)
    # Output theoretical tip deflection
    print(tipDeflection)
    
    # Deform the beam
    err, res = testSprg.forward_integration(testSprg.deform_ODE, np.array([0,30,0]), 0)

    # Output simulator tip deflection in y direction
    print(res[2,-1])

    assert np.isclose(res[2,-1], tipDeflection, rtol=0.01) # relatively large tolerance due to hard assumptions in EB beam theory 

    # print("In this plot you should see a tip deflection of 1 inch")
    # testSprg.plot_deform(testSprg)

# @pytest.mark.skip()
def test_moment_uniform_curved_beam():

    """ This test ensures that a curved beam with a profile simple enough to yield an algebraic solution
        agrees with theory when bent by a pure moment under the simulator """


    # Numerically determined parameters to result in unit displacement
    R = 1.4564
    M = 67.81002
    Ic = 0.005
    # Define a spring with constant curvature and constant Ic (therefore constant thickness)
    testPath = LinearRnSpiral(1, R*np.pi, startPoint=(R,0), endPoint=(-R,0),
                                       initialRadius=R,  finalRadius=R)
    testCrsc = Constant_Ic(testPath, 0.375, Ic0 = Ic)
    testSprg = Spring(testPath, testCrsc, TestMaterial, resolution = 200)

    # Determine predicted deflection using curved beam theory
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

    # Use simulator to deflect beam
    err, res = testSprg.forward_integration(testSprg.deform_ODE, np.array([0,0,M]), M)

    DX = res[1,-1]-testSprg.xL
    DY = res[2,-1]-testSprg.yL

    print("Numerically calculated displacement")
    print(DX, DY)
    print(testSprg.dist)

    # testSprg.plot_deform(False)

    assert(np.isclose(testSprg.dist, predictDist))
    assert(np.isclose(DX, predictDX))
    
# @pytest.mark.skip()
def test_Yforce_uniform_curved_beam():
    
    """ This test is meant to ensure the beam simulator agrees with theory for a curved beam under tip force load. 
        For tip force loads, curved beams do not yield closed-form solutions for bending, and therefore this test simoly
        solves a simplified version of the bending problem numerically. Basically, this test confirms that the general 
        solver reduces to a reasonable solution in the case of a tip force loading in the Y direction. """

    # Define a spring with constant curvature and constant Ic (therefore constant thickness)
    R = 1
    Ic = 0.01
    
    testPath = LinearRnSpiral(1, R*np.pi, startPoint=(R,0), endPoint=(-R,0),
                                       initialRadius=R,  finalRadius=R)
    testCrsc = Constant_Ic(testPath, 0.375, Ic0 = Ic)
    testSprg = Spring(testPath, testCrsc, TestMaterial, resolution = 200)

    # Apply force in the Y direction
    Fy = -TestMaterial.E*Ic/10
    Fx = 0
    
    # Use simulator to deform beam
    err, res = testSprg.forward_integration(testSprg.deform_ODE, np.array([0,Fy,0]), 0)

    # There is no closed form algebraic solution for gamma, x, y, so we instead integrate a simplified form of the general ODE for this case
    def simplified_ODE(coord, states):
        s = coord
        gamma, x, y = states
        LHS = np.empty(3)
        LHS[0] = -Fy/(TestMaterial.E*Ic)*(x+R)
        LHS[1] = np.cos(s/R+testPath.get_alpha(0)+gamma)
        LHS[2] = np.sin(s/R+testPath.get_alpha(0)+gamma)

        return LHS

    # Use simplified ODE to deform beam
    SimplifiedSolution = utils.fixed_rk4(simplified_ODE, np.array([0,R,0]), testSprg.ximesh)
    GSimplified = SimplifiedSolution[0,:]
    XSimplified = SimplifiedSolution[1,:]
    YSimplified = SimplifiedSolution[2,:]
    
    print("")
    print("Attained Result:", res[:,-1])
    print("Expected Result:", SimplifiedSolution[:,-1])

    assert(np.all(np.isclose(res[:,-1], SimplifiedSolution[:,-1])))
    
    # testSprg.plot_deform(False)

def test_Xforce_uniform_curved_beam():
    
    """ This test is meant to ensure the beam simulator agrees with theory for a curved beam under tip force load. 
        For tip force loads, curved beams do not yield closed-form solutions for bending, and therefore this test simoly
        solves a simplified version of the bending problem numerically. Basically, this test confirms that the general 
        solver reduces to a reasonable solution in the case of a tip force loading in the X direction. """

    # Define a spring with constant curvature and constant Ic (therefore constant thickness)
    R = 1
    Ic = 0.01
    
    testPath = LinearRnSpiral(1, R*np.pi, startPoint=(R,0), endPoint=(-R,0),
                                       initialRadius=R,  finalRadius=R)
    testCrsc = Constant_Ic(testPath, 0.375, Ic0 = Ic)
    testSprg = Spring(testPath, testCrsc, TestMaterial, resolution = 200)

    # Apply force in the X direction
    Fy = 0
    Fx = -TestMaterial.E*Ic/10

    # Use simulator to deform beam
    err, res = testSprg.forward_integration(testSprg.deform_ODE, np.array([Fx,0,0]), 0)

    # There is no closed form algebraic solution for gamma, x, y, so we instead integrate a simplified form of the general ODE for this case
    def simplified_ODE(coord, states):
        s = coord
        gamma, x, y = states
        LHS = np.empty(3)
        LHS[0] = Fx/(TestMaterial.E*Ic)*(y)
        LHS[1] = np.cos(s/R+testPath.get_alpha(0)+gamma)
        LHS[2] = np.sin(s/R+testPath.get_alpha(0)+gamma)

        return LHS

    # Use simplified ODE to deform beam
    SimplifiedSolution = utils.fixed_rk4(simplified_ODE, np.array([0,R,0]), testSprg.ximesh)
    GSimplified = SimplifiedSolution[0,:]
    XSimplified = SimplifiedSolution[1,:]
    YSimplified = SimplifiedSolution[2,:]
    
    print("")
    print(testSprg.dgds0)
    print("Attained Result:", res[:,-1])
    print("Expected Result:", SimplifiedSolution[:,-1])
    assert(np.all(np.isclose(res[:,-1], SimplifiedSolution[:,-1])))