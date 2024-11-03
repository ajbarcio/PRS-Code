import numpy as np
import json

import os
import subprocess
import argparse

import manualSpringDesigner
from manualSpringDesigner import defineSpring

from matplotlib import pyplot as plt

testingSpring = defineSpring()

torques = np.linspace(-testingSpring.torqueCapacity, testingSpring.torqueCapacity, 16)
deflections = []
for torque in torques:
    res, SF, divergeFlag, i = testingSpring.deformMode(torque,testingSpring.deform_ODE)
    deflections.append(testingSpring.dBeta)

plt.plot(torques, deflections)
plt.show()