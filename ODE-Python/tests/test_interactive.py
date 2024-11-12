import pytest
import inspect

import numpy as np
from numpy import linalg as lin
from matplotlib import pyplot as plt
import filecmp
import random
from unittest.mock import MagicMock
import os

from modules.PATHDEF import Path, TestPath, LinearRnSpiral, RadiallyEndedPolynomial
from modules.CRSCDEF import Crsc, Constant_Ic, Piecewise_Ic_Control
from modules.materials import TestMaterial
from modules.spring import Spring, determineFastestSolver
from modules.interactive import Interactive_Spring

from modules.utils import deg2rad

from modules.StatProfiler import SSProfile

mock_event = MagicMock()

# @pytest.mark.skip(reason="Interactive deprecated")
def test_interactive():
    print("\n")
    # Initialize objects for test
    testPath = RadiallyEndedPolynomial(2,6)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([.008,.0008,.008]))
    testSprg = Spring(testPath, testCrsc, TestMaterial)
    # Initialize Interactive Spring for test
    testInteractive = Interactive_Spring(testPath, testCrsc, TestMaterial, name="Test")

    # Test basic "sliders change intrinsic spring values" functionality
    # THIS ONLY WORKS FOR THE POLYNOMIAL SPRING, CURRENTLY
        # This also tests the update plot callback

    # Get some initial param and surface files
    params1 = testInteractive.export_parameters()
    surfcs1 = testInteractive.export_surfaces(mock_event)
    
    # Screw with the values in a random way
    for slider in testInteractive.adjusters:
        slider.set_val(random.random()*(slider.valmax-slider.valmin)+slider.valmin)

    # Get the final param and surface files
    params2 = testInteractive.export_parameters()
    surfcs2 = testInteractive.export_surfaces(mock_event)

    # Make sure they changed
    for a, b in zip(params1, params2):
        assert not filecmp.cmp(a, b, shallow=False)
    for a, b in zip(surfcs1, surfcs2):
        assert not filecmp.cmp(a, b, shallow=False)

    for a in params1+params2+surfcs1+surfcs2:
        try:
            os.remove(a)
            # print(f"File {a} removed successfully.")
        except FileNotFoundError:
            print(f"File {a} not found.")
        except PermissionError:
            print(f"Permission denied to remove {a}.")
        except Exception as e:
            print(f"Error occurred: {e}")

    return True

@pytest.mark.skip(reason="This functionality is in development")
def test_loadValues():
    print("\n")
    # Initialize objects for test
    testPath = RadiallyEndedPolynomial(2,6)
    testCrsc = Piecewise_Ic_Control(testPath, IcPts=np.array([.008,.0008,.008]))
    testSprg = Spring(testPath, testCrsc, TestMaterial)
    # Initialize Interactive Spring for test
    testInteractive = Interactive_Spring(testPath, testCrsc, TestMaterial, name="Test")

    # Test basic "sliders change intrinsic spring values" functionality
    # THIS ONLY WORKS FOR THE POLYNOMIAL SPRING, CURRENTLY
        # This also tests the update plot callback

    # Get some initial param and surface files
    params1 = testInteractive.export_parameters()
    surfcs1 = testInteractive.export_surfaces(mock_event)
    
    # Screw with the values in a random way
    for slider in testInteractive.adjusters:
        slider.set_val(random.random()*(slider.valmax-slider.valmin)+slider.valmin)

    # Get the final param and surface files
    params2 = testInteractive.export_parameters()
    surfcs2 = testInteractive.export_surfaces(mock_event)

    # Make sure they changed
    for a, b in zip(params1, params2):
        assert not filecmp.cmp(a, b, shallow=False)
    for a, b in zip(surfcs1, surfcs2):
        assert not filecmp.cmp(a, b, shallow=False)

    # for 
    # testSprg.load_values()

    for a in params1+params2+surfcs1+surfcs2:
        try:
            os.remove(a)
            # print(f"File {a} removed successfully.")
        except FileNotFoundError:
            print(f"File {a} not found.")
        except PermissionError:
            print(f"Permission denied to remove {a}.")
        except Exception as e:
            print(f"Error occurred: {e}")

    return True