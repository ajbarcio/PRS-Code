import numpy as np
import numpy.linalg as lin

# err = 0.25
# J = np.array([0,1,2,3])
# Jinv = lin.pinv(np.atleast_2d(J).T)
# print(Jinv)
# print(Jinv*err)
# print(J+Jinv)

Aname = "A_surface.txt"
Bname = "B_surface.txt"

A = np.genfromtxt("surfaces\\"+Aname, delimiter=',')
B = np.genfromtxt("surfaces\\"+Bname, delimiter=',')
springName = "manual_spring"

A = A[:,0:2]
A = list(map(tuple,A))
print(A)

params = np.genfromtxt("springs\\"+springName, delimiter=',')
print(params)