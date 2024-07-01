import numpy as np
from scipy.integrate import solve_ivp as variable_mesh_solve
import numpy.linalg as lin
import matplotlib.pyplot as plt
import time
from spring import Spring
import spring
from materials import Maraging300Steel
from StatProfiler import SSProfile
from utils import fixed_rk4

deg2rad = np.pi/180
# Give the spring some parameters that hopefully do not affect anything about the math
defaultSpring = spring.generate_simple_spring()
defaultSpring.full_results(deformBool=False)
# Rootfinding method for getting outer surface geometry
defaultSpring.generate_undeformed_surfaces()
print(defaultSpring.rootfindingSteps)
# Solve the governing equations to get the implied input function values
# This method is a bit incomplete, as it assumes rn agrees where used as an input to the method
# However, since it provides low error, it is probably fine
checkIc = defaultSpring.t*defaultSpring.rn*(defaultSpring.lb+defaultSpring.la)*(defaultSpring.lb/2-defaultSpring.la/2)
checkrn = (defaultSpring.la+defaultSpring.lb)/(np.log((defaultSpring.rn+defaultSpring.lb)/(defaultSpring.rn-defaultSpring.la)))
# Check that the input and resultant rn, cI agree
plt.figure("Checking rootfinding accuracy")
plt.plot(defaultSpring.smesh, checkIc-defaultSpring.Ic, label="Ic error")
plt.plot(defaultSpring.smesh, checkrn-defaultSpring.rn, label="rn error")
plt.figure("raw")
plt.plot(defaultSpring.smesh, defaultSpring.rn)
plt.plot(defaultSpring.smesh, defaultSpring.Ic)
plt.figure("la, lb")
plt.plot(defaultSpring.smesh, defaultSpring.la)
plt.plot(defaultSpring.smesh, defaultSpring.lb)
plt.legend()

la_lb_0 = np.array([defaultSpring.la[50], defaultSpring.lb[50]])
result  = fixed_rk4(defaultSpring.geo_ODE, la_lb_0, defaultSpring.ximesh[50:])
print(result)
plt.plot(defaultSpring.ximesh[50:], result[0])
plt.plot(defaultSpring.ximesh[50:], result[1])
plt.show()