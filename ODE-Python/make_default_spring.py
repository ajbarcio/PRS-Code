import numpy as np
import numpy.linalg as lin
import matplotlib.pyplot as plt
from os import remove
import glob

from spring import Spring
import spring
from materials import Maraging300Steel

defaultSpring = spring.generate_default_spring()
defaultSpring.full_results(deformBool=False)

A = np.hstack((defaultSpring.undeformedASurface,np.atleast_2d(np.zeros(len(defaultSpring.undeformedASurface))).T))
B = np.hstack((defaultSpring.undeformedBSurface,np.atleast_2d(np.zeros(len(defaultSpring.undeformedBSurface))).T))
np.savetxt("surfaces\\default_A_surface.txt", A, delimiter=",", fmt='%f')
np.savetxt("surfaces\\default_B_surface.txt", B, delimiter=",", fmt='%f')

plt.show()
