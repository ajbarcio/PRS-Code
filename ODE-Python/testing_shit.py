import numpy as np
import scipy as sp
import numpy.linalg as lin
import scipy.linalg as slin
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *
from openpyxl import *



p6 = np.array([.5, 2])
p3 = np.array([1, 1])
p9 = np.array([-2, -1])

#print(4*p6)
A = np.array([[-1,4,0,1,0,0],[0,1,1,0,0,0],[0,0,-1,4,0,1],[0,0,0,1,1,0]])
B = np.array([4*p3,2*p3,4*p6,2*p6])
print(A)
print(B)
#print(A.dot(A.T))
#print(lin.det(A.dot(A.T)))
print(A.T.dot((lin.inv(A.dot(A.T))).dot(B)))
fun = lambda X: A.dot(X)-B

con1fun = lambda X: X[0,0]
con1 = sp.optimize.NonlinearConstraint(con1fun, 0, 0)
con2fun = lambda X: np.tan(X[5, 1]/X[5, 0])
con2 = sp.optimize.NonlinearConstraint(con2fun, np.tan(p9[1]/p9[0]), np.tan(p9[1]/p9[0]))

cons = (con1, con2)

X0 = A.T.dot((lin.inv(A.dot(A.T))).dot(B))

result = sp.optimize.minimize(fun, X0, constraints=cons)



#result = lin.solve(A.dot(A.T),B.T.dot(A))
#print(result)

