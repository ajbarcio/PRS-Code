import numpy as np

deg2rad = np.pi/180

Rin = 1
Rout = 6
dBMax = 5*deg2rad
thmin = 30*deg2rad
thmax = 360*deg2rad
q = 1.4
c = (Rout-Rin)/thmax**q

def spiral_radius(theta, R):
    R = c*theta^q + Rin

