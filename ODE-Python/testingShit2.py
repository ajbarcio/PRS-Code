from stiffness_library import *

deg2rad = np.pi/180
n = 2
fullArcLength = 5

E = 27.5*10**6
outPlaneThickness = .375

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2

R1 = (R0+R3)/2
R2 = (R0+R3)/2

betaB = 160/3*deg2rad*.5
betaC = 2*160/3*deg2rad
betaD = 160*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.008, .001, .018])

ctrlcIs      = np.array([0,fullArcLength*.5,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

# maxTorque = 4554.5938
maxTorque = 5000
maxDBeta  = 0.087266463

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]

geometryDef , smesh = drag_vector_spring(dragVector0)
la,lb = outer_geometry(5,geometryDef,False)
print(la,lb)