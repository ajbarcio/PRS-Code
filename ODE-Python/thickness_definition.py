import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from StatProfiler import SSProfile
import scipy
from copy import deepcopy as dc

deg2rad = np.pi/180

class Piecewise_Ic_Control():
    def __init__(self, 