import numpy as np
import numpy.linalg as lin

from matplotlib import pyplot as plt

from CRSCDEF import Constant_Ic
from PATHDEF import LinearRnSpiral, RadiallyEndedPolynomial
from materials import TestMaterial
from utils import PPoly_Eval
from spring import Spring

# path = LinearRnSpiral(2, None, startPoint=(1,1), endPoint=(-1,3))
path = RadiallyEndedPolynomial(2, 6)
crsc = Constant_Ic(path, 0.5, Ic0=0.0125)
mtrl = TestMaterial
sprg = Spring(path, crsc, TestMaterial)
sprg.plot_spring()
la, lb = crsc.get_neutralDistances(sprg.resl)
plt.show()