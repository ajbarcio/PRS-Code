import numpy as np
from spring2 import Spring
import materials
import PATHDEF
import CRSCDEF

# Standard procedure: define a path, then define a thickness profile using the
#                     length of that path, then pass both of those objects into
#                     the overall spring object that contains deformation 
#                     methods

defaultPath = PATHDEF.Minimal_Polynomial_Definition()
defaultGeom = CRSCDEF.Piecewise_Ic_Control(path=defaultPath)

defaultSpring = Spring(defaultGeom, materials.Maraging300Steel())

