Ic = 0.002
rn = 1


# treating path as known (nominal)
# treating IC   as known (preliminary)
Ic = self.crsc.get_Ic(xi)
rn = self.path.get_rn(xi)

dIcdxi = self.crsc.get_dIc(xi)
drndxi = self.path.get_drn(xi)

# # transform to arc length space

# dxids = self.path.get_dxi_n(xi)
# dIcds = dIcdxi*dxids
# drnds = drndxi*dxids

# Version without transform

dIcds = dIcdxi
drnds = drndxi

# q = prime     states (la, lb)
# p = secondary states (Ic, rn)
# Q = assembled prime     quasi-state vectors
# P = assembled secondary quasi-state vectors

Q = np.array([[-q[0],            q[1]],
            [1/(rn-q[0])-1/rn, 1/(rn+q[1])-1/rn]])
P = np.array([[1/(rn*self.t),    -Ic/(rn**2*self.t)],
            [0,                1/(rn-q[0])-1/(rn+q[1])-(q[1]+q[0])/rn**2]])
p_dot = np.array([[dIcds],[drnds]])

# This is the diffeq step

try:
    q_dot = lin.inv(Q).dot(P).dot(p_dot)
except lin.LinAlgError:
    print("excepted", xi)
    q_dot = lin.pinv(Q).dot(P).dot(p_dot)

# Back to xi space
# q_dot = q_dot/dxids

return q_dot.flatten()