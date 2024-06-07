import numpy as np
import numpy.linalg as lin

err = 0.25
J = np.array([0,1,2,3])
Jinv = lin.pinv(np.atleast_2d(J).T)
print(Jinv)
print(Jinv*err)
print(J+Jinv)

import warnings

def fxn():
    warnings.warn("fuck you", DeprecationWarning)
    print("called")
    
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    fxn()
