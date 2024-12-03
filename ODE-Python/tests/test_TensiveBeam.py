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


def test_tensive_straight_beam():
    testPath = RadiallyEndedPolynomial(1, 10, radii=np.array([0, 10]), ffradii=np.array([0, 10]), alphaAngles=np.array([0,0]), betaAngles=np.array([0,0]), XYFactors=np.array([]))
    # testCrsc = Constant_Ic(testPath, 0.375, Ic0=0.00125)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([1/120,1/120]), t=0.1)
    testSprg = Spring(testPath, testCrsc, TestMaterial)
    print("hey!")
    print(testPath.get_rn(testSprg.ximesh))

    testSprg.forward_integration(testSprg.deform_withTension_ODE, np.array([1000,0,0]), 0)
    testSprg.plot_deform(testSprg)

def main():
    test_tensive_straight_beam()

if __name__ == "__main__":
    main()