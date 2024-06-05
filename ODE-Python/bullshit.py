import numpy as np
import numpy.linalg as lin

err = 0.25
J = np.array([0,1,2,3])
Jinv = lin.pinv(np.atleast_2d(J).T)
print(Jinv)
print(Jinv*err)
print(J+Jinv)