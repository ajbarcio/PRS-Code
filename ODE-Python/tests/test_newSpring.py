import pytest
import inspect

import numpy as np
from numpy import linalg as lin
from matplotlib import pyplot as plt
import filecmp
import random
from unittest.mock import MagicMock
import os

from modules.PATHDEF import Mirabilis
from modules.CRSCDEF import Constant_Ic
from modules.materials import TestMaterial
from modules.spring import Spring, determineFastestSolver, Spring2
from modules.interactive import Interactive_Spring

from modules.utils import deg2rad

from modules.StatProfiler import SSProfile

def test_forwardIntegrationWithSciPy():
    path = Mirabilis(2, 1.1, 2.5, np.pi)
    crsc = Constant_Ic(path, 0.375, Ic0=0.01318)
    spring = Spring2(path, crsc, TestMaterial, name="newSpringTest")
    arcLens, states, sol = spring.spring_forward_solve(spring.ODE_S, [50,50,50])
    smesh = np.linspace(0,spring.path.arcLen,100)
    plt.plot(smesh, sol(smesh).T)

def test_plotting():
    path = Mirabilis(2, 1.1, 2.5, np.pi)
    crsc = Constant_Ic(path, 0.375, Ic0=0.01318)
    spring = Spring2(path, crsc, TestMaterial, name="newSpringTest")
    smesh = np.linspace(0,spring.path.arcLen,100)
    ximesh = path.theta(smesh)
    print(crsc.get_Thk(ximesh))
    spring.plot_spring(100, showBool=True)