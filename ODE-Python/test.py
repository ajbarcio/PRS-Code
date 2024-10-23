import numpy as np
from numpy import linalg as lin
from scipy import optimize as opt
from matplotlib import pyplot as plt
import scipy.special as spec
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

from scipy.optimize import fsolve
designStress = 180000
rn0 = 5
E = 27500000
M = 2500
Fx = 0
Fy = 0
n=2
t=0.375
momentArmY=2
momentArmX=1

la0 = np.linspace(0,0.5,100)
lb0 = lb0   = rn0*-spec.lambertw((np.exp(la0/rn0-1)*(la0-rn0))/rn0,-1)-rn0
lb0 = np.real(lb0)
Ic0   = t*rn0/2*(lb0**2-la0**2)
dgds0 = (M/n + momentArmY*Fx - momentArmX*Fy)/(E*Ic0)

designStress = E*rn0*dgds0*np.maximum(la0/(rn0-la0),lb0/(rn0+lb0))

plt.figure("solution")
plt.plot(la0,lb0)
plt.figure("difference")
plt.plot(la0,lb0-la0)
plt.figure("root stress")
plt.plot(la0,designStress)
plt.plot(la0,180000*np.ones_like(la0))

def findTargetRootStress(w):
    la = w
    F = np.array(1)
    lb = rn0*-spec.lambertw((np.exp(la/rn0-1)*(la-rn0))/rn0,-1)-rn0
    lb = np.real(lb)
    Ic = self.t*rn0/2*(lb**2-la**2)
    dgds0 = 
    F = self.E*rn0*dgds0*max(la0/(rn0-la),lb/(rn0+lb))-self.designStress
    return F
la0   = .32
solutionInfo = fsolve(findTargetRootStress,la0,full_output=1)
solution = solutionInfo[0]
la0 = solution[0]
lb0   = rn0*-spec.lambertw((np.exp(la0/rn0-1)*(la0-rn0))/rn0,-1)-rn0
lb0 = np.real(lb0)
print(lb0)


plt.show()





# def nonlinear(w):
#     la, lb, dgds0, Ic0 = w
#     F=np.zeros(4)
#     F[0] = designStress-E*rn*dgds0*la/(rn-la)
#     F[1] = rn-(la+lb)/np.log((rn+lb)/(rn-la))
#     F[2] = dgds0-(M/n+Fx*momentArmY-Fy*momentArmX)/(E*Ic0)
#     F[3] = Ic0-t*rn/2*(lb**2-la**2)
#     return F
# # generate an initial guess
# initialGuess=np.random.rand(4)    
 
# # solve the problem    
# solutionInfo=fsolve(nonlinear,initialGuess,full_output=1)
# print(solutionInfo)





# W = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
# Q = np.array([[1,2],[3,4]])
# P = np.array([[1,0],[-1,1]])
# qStar = np.array([[1],[2]])
# pStar = np.array([[3],[0]])

# QP = np.vstack((Q,P))
# qStarPStar = np.vstack((qStar,pStar))
# qdot = lin.pinv(W.dot(QP)).dot(W).dot(qStarPStar)
# print(qdot)
# print(qdot[0])
# print(qdot[1])

# class TestClass:
#     def __init__(self):
#         self.a = 2
#         self.b = 3

# test = TestClass()
# test.define_some_stuff()
# # test.define_some_stuff.add_stuff()
# print(test.c)

# x = [0, 1]
# y = [0, 1]
# X, Y = np.meshgrid(x, y)
# z = np.zeros_like(X)
# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})

# ax.plot_surface(X,Y,z)
# plt.show()

# states = [0,1,2,3,4]
# gamma, x, y, la, lb = states
# print(gamma)
# print(la)

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
# X = X.subs(lb,la)
# print(X[0])
# print(X[1])

# list  = [0,3,5,4,8,6,1,2,3,2,7,9]
# array = np.array(list)
# elist = enumerate(list)
# earray = np.ndenumerate(array)
# # print(earray[:])
# limit = 3
# filteredList  = [item for item in list if item < limit]
# filteredArray = [item for item in array if item < limit]

# filteredEList  = [pair for pair in elist if pair[1] < limit]
# filteredEArray = [pair for pair in earray if pair[1] < limit]

# print(filteredList)
# print(filteredArray)

# print(filteredEList)
# print(filteredEArray)
# dt = np.dtype('float,float')
# print(np.array([i[0] for i in filteredEArray]))

# a = [0,1,2]
# b = [4,1,4]
# c = [2,3,4]

# bigArray = np.vstack((a,b,c))
# print(bigArray)

# filteredBigArray = [column for column in bigArray.T if any(column[1:] > limit)]
# print(filteredBigArray)

# repairedArray = np.array(filteredBigArray).T
# print(repairedArray)

# repairedIndices = repairedArray[0,:]
# print(repairedIndices)

# d = [5,4,3,2,1]
# print([d[i] for i in repairedIndices])
# print([])