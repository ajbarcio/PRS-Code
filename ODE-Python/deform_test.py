import numpy as np

import materials
import PATHDEF
import CRSCDEF
from spring import Spring, deg2rad

straightPath = PATHDEF.Minimal_Polynomial_Definition4(n=1, fullParamLength=1,
                                                  radii = np.array([1,2]),
                                                  ffradii = np.array([1,2]),
                                                  alphaAngles = np.array([0,0]),
                                                  betaAngles = np.array([0,0]),
                                                  XYFactors=np.array([]))

Ic = 1/12*.375*.125**3

straightCrsc = CRSCDEF.Piecewise_Ic_Control(pathDef=straightCrsc,
                                        IcPts = np.array([Ic, Ic]),
                                        IcParamLens=np.array([]))
straightPath.get_crscRef()

straightSpring = Spring(straightCrsc, materials.Maraging300Steel,
                        name="20270730_spring")

