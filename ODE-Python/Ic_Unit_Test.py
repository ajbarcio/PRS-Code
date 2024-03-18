import stiffness_library as st
import numpy as np
import matplotlib.pyplot as plt

deg2rad = np.pi/180

globalInnerRadiusLimit = 0.75
globalOuterRadiusLimit = 6/2

R0 = 2.2/2
R3 = 5.9/2*.9

R1 = (R0+R3)/2+.26
R2 = (R0+R3)/2-.25

betaB = 150/3*deg2rad*.5
betaC = 2*150/3*deg2rad
betaD = 150*deg2rad
beta0 = betaD

x0 = R0
y0 = 0

fullArcLength = 5.2

pts = np.array([[x0, y0],[R1*np.cos(betaB),R1*np.sin(betaB)],[R2*np.cos(betaC),R2*np.sin(betaC)],[R3*np.cos(betaD),R3*np.sin(betaD)]])
cIs  = np.array([.008, .001, .008])

ctrlcIs      = np.array([0,fullArcLength*.5,fullArcLength])
ctrlLengths = np.array([0,fullArcLength*0.333,fullArcLength*0.667,fullArcLength])

maxTorque = 13541.64
maxDBeta  = 0.087266463

dragVector0 = [R0, R1, R2, R3, \
               betaB, betaC, beta0, \
               cIs[0], cIs[1], cIs[2], \
               ctrlcIs[1], \
               ctrlLengths[1], ctrlLengths[2],
               fullArcLength]

cICoeffs = st.cI_poly(cIs, ctrlcIs)
smesh = np.linspace(0,fullArcLength,100)
assert(st.cI_s(0,cICoeffs)==cIs[0])
print(st.cI_s(0,cICoeffs))
assert(round(st.cI_s(ctrlcIs[1], cICoeffs), 3)==cIs[1])
print(round(st.cI_s(ctrlcIs[1], cICoeffs), 3))
assert(round(st.cI_s(fullArcLength,cICoeffs),3)==cIs[2])
print(st.cI_s(fullArcLength,cICoeffs))
plt.plot(st.cI_s(smesh, cICoeffs))
plt.show()