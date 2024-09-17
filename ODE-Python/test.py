import numpy as np
from numpy import linalg as lin
from scipy import optimize as opt
from matplotlib import pyplot as plt

import pandas as pd

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

from utils import fixed_rk4, numerical_fixed_mesh_diff

import sympy as sp

df = pd.read_excel('Spring_Constraints.ods', engine='odf', index_col=0)

data =df.loc['Size 5', 'ID lim (in)']
print(data)