import numpy as np
from numpy import linalg as lin
from scipy import optimize as opt
from matplotlib import pyplot as plt
import os
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
# x,y = sp.symbols("x,y")
# la, lb, rn, drn = sp.symbols("la,lb,rn,drn")
# aleph = 1/(rn+x)-1/rn
# beth  = 1/(rn-x)-1/(rn+y)-(x+y)/rn**2

# # print(aleph.subs(x,la))
# # print(beth.subs([(x,y),(la,lb)]))

# A = sp.Matrix(2,2,[-la,lb,aleph.subs(x,-la),aleph.subs(x,lb)])
# B = sp.Matrix(2,1,[0,drn*beth.subs([(x,la),(y,lb)])])
# X = A.inv()*B
# print(X[0])
# print(X[1])

list  = [0,3,5,4,8,6,1,2,3,2,7,9]
array = np.array(list)
elist = enumerate(list)
earray = np.ndenumerate(array)

for a, b in enumerate(list):
    print(a,b)
# print(earray[:])
limit = 3
filteredList  = [item for item in list if item < limit]
filteredArray = [item for item in array if item < limit]

filteredEList  = [pair for pair in elist if pair[1] < limit]
filteredEArray = [pair for pair in earray if pair[1] < limit]

print(filteredList)
print(filteredArray)

print(filteredEList)
print(filteredEArray)
dt = np.dtype('float,float')
print(np.array([i[0] for i in filteredEArray]))

a = [0,1,2]
b = [4,1,4]
c = [2,3,4]

bigArray = np.vstack((a,b,c))
print(bigArray)

filteredBigArray = [column for column in bigArray.T if any(column[1:] > limit)]
print(filteredBigArray)

repairedArray = np.array(filteredBigArray).T
print(repairedArray)

repairedIndices = repairedArray[0,:]
print(repairedIndices)

d = [5,4,3,2,1]
print([d[i] for i in repairedIndices])