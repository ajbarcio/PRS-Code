import numpy as np
from spring import Spring
import materials
import PATHDEF
import CRSCDEF

# Standard procedure: define a path, then define a thickness profile using the
#                     length of that path, then pass both of those objects into
#                     the overall spring object that contains deformation 
#                     methods

testPath = PATHDEF.Minimal_Polynomial_Definition()
testCrsc = CRSCDEF.Piecewise_Ic_Control(path=testPath)
# fuck this
testPath.get_crscRef(testCrsc)


defaultSpring = Spring(defaultGeom, materials.Maraging300Steel())

