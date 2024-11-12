import pytest
import inspect

import numpy as np
from numpy import linalg as lin
from matplotlib import pyplot as plt
import filecmp
import random
from unittest.mock import MagicMock, Mock
import os
from scipy.integrate import solve_ivp

from modules.PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from modules.CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control
from modules.materials import TestMaterial
from modules.spring import Optimized_Spring, determineFastestSolver
from modules.interactive import Interactive_Spring

from modules.utils import deg2rad

from modules.StatProfiler import SSProfile

@pytest.mark.skip(reason="takes forever and isnt strictly necessary")
def test_weight_tuning():
    # Initialize objects for test
    testPath = RadiallyEndedPolynomial(2,6)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([.008,.0008,.008]))
    testSprg = Optimized_Spring(testPath, TestMaterial)

    threshold, transition = testSprg.tune_weightings()
    print(threshold, transition)

# @pytest.mark.skip()
def test_optimization_class():
    # Initialize objects for test
    testPath = RadiallyEndedPolynomial(2,6)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([.008,.0008,.008]))
    testSprg = Optimized_Spring(testPath, TestMaterial)

    # Make sure this function runs without errors
    testSprg.prepare_weightings()

    # test forward integration in some small pure moment case
    # err, res = testSprg.forward_integration(testSprg.weighted_ODE, np.array([-50,-50,testSprg.torqueCapacity]), 0)
    h0 = np.sqrt(6*testSprg.torqueCapacity/testSprg.n/(testSprg.t*testSprg.designStress))
    Ic0 = testSprg.t*(h0)**3/12
    print(h0)
    la0=h0/2
    lb0=la0
    # print("init thickness:", la0)
    print(Ic0)
    res = solve_ivp(testSprg.weighted_ODE, 
                    [0, testSprg.fullArcLength], 
                    [0,testSprg.x0,testSprg.y0,la0,lb0], 
                    args=[0,0,testSprg.torqueCapacity])
    
    plt.plot(res.t, res.y.T[:,0])
    plt.figure()
    plt.plot(res.t, res.y.T[:,1])
    plt.plot(res.t, res.y.T[:,2])
    plt.figure()
    plt.plot(res.t, res.y.T[:,3])
    plt.plot(res.t, res.y.T[:,4])
    plt.show()

    # plt.plot(testSprg.xiData, testSprg.stressData, label="stress")
    # plt.plot(testSprg.xiData, [i*testSprg.designStress for i in testSprg.weightData], label="curved behavior weight")
    # plt.plot(testSprg.xiData, [(i+1-i)*testSprg.designStress for i in testSprg.weightData], label = "design stress")
    # plt.legend()
    # plt.figure()
    
    # la = res[3,:]
    # lb = res[4,:]
    # h = la+lb

    # plt.plot(testSprg.ximesh, la)
    # plt.plot(testSprg.ximesh, lb)
    # plt.plot(testSprg.ximesh, h)
    
    # rn = testSprg.path.get_rn(testSprg.ximesh)
    # print(rn[0])
    # Ic = testSprg.t*rn/2*(lb**2-la**2)
    
    # plt.figure()
    # plt.plot(testSprg.ximesh, Ic)
    
    
    # plt.show()
    # # # test forward integration in some reasonable forced case
    # # err, res = testSprg.forward_integration(testSprg.weighted_ODE, np.array([testSprg.torqueCapacity/20,-testSprg.torqueCapacity/60,testSprg.torqueCapacity]), 0)


    # # Test degenerate case
    # testSprg.deform_by_torque(0, testSprg.weighted_ODE,
    #                           np.array([0,0,0]))
