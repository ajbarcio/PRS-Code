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

# df = pd.read_excel('Spring_Constraints.ods', engine='odf', index_col=0)

# data =df.loc['Size 5', 'ID lim (in)']
# print(data)

# a,b,c,d,e = sp.symbols("a,b,c,d,e")
x,y = sp.symbols("x,y")
la, lb, rn, drn = sp.symbols("la,lb,rn,drn")
aleph = 1/(rn+x)-1/rn
beth  = 1/(rn-x)-1/(rn+y)-(x+y)/rn**2

# print(aleph.subs(x,la))
# print(beth.subs([(x,y),(la,lb)]))

A = sp.Matrix(2,2,[-la,lb,aleph.subs(x,-la),aleph.subs(x,lb)])
B = sp.Matrix(2,1,[0,drn*beth.subs([(x,la),(y,lb)])])
X = A.inv()*B
print(X[0])
print(X[1])