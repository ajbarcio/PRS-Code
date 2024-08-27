import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt

from materials import Maraging300Steel
from utils import fixed_rk4, numerical_fixed_mesh_diff, colorline

rn = 1.2
Ic = .005
# get derivatives in terms of parameter
dIcdxi = -.00005
drndxi = -.5
# convert to arc length space
dxids = 0.99

dIcds = dIcdxi*dxids
drnds = drndxi*dxids

t = .375

q = np.array([0.125,0.123])

# This is some good good math I swear
P = np.array([[1/(rn*t),    -Ic/(rn**2*t)],
                [0,         1/(rn-q[0])-1/(rn+q[1])-(q[0]+q[1])/rn**2]])
print(P)

Q = np.array([[-q[0],            q[1]],
              [1/(rn-q[0])-1/rn, 1/(rn+q[1])-1/rn]])
print(Q)
p_dot = np.array([[dIcds],[drnds]])
print(p_dot)
q_dot = lin.inv(Q).dot(P).dot(p_dot)
print(q_dot)