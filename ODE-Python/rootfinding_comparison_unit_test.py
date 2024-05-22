import numpy as np
from scipy.integrate import solve_ivp as variable_mesh_solve
import numpy.linalg as lin
import matplotlib.pyplot as plt
import time
from stiffness_library import *
from StatProfiler import SSProfile

# Give the spring some parameters that hopefully do not affect anything about the math


deg2rad = np.pi/180
# dragVector0 = [1, 2, 2.5, 4, \
#                5*deg2rad, 10*deg2rad, 15*deg2rad,  \
#                .001, .0005, .001, \
#                fullArcLength/2,  fullArcLength/3, 2*fullArcLength/3, fullArcLength]
dragVector0 = [1.10000000e+00, 2.13750000e+00, 1.62750000e+00, 2.65500000e+00,
 4.36332313e-01, 1.74532925e+00, 2.79252680e+00, 4.74085857e-03,
 7.40837797e-04, 4.74084339e-03, 2.60000000e+00, 1.73160000e+00,
 3.46840000e+00, 5.20000000e+00]
fullArcLength = dragVector0[-1]
# dragVector0 = [1.0, 2.0, 6.0, 10.0, 1*deg2rad, 2*deg2rad, 3*deg2rad, 0.5, 0.5, 0.5, \
            #    6, 4, 8, 12]
print("Spring Parameters:", dragVector0)

# Create the spring geometry
geometryDef, smesh = drag_vector_spring(dragVector0)
# Create a mesh size
smesh = np.linspace(0, fullArcLength, 1001)
# Get vectors for our "input" functions
rn = np.empty(len(smesh))
Ic = np.empty(len(smesh))
for i in range(len(smesh)):
    rn[i] = r_n(smesh[i], geometryDef[0], geometryDef[1])
    Ic[i] = PPoly_Eval(smesh[i], geometryDef[2])

# Rootfinding method for getting outer surface geometry
print("start lAB method")
start = time.time()
la = np.empty(len(smesh))
lb = np.empty(len(smesh))
hLALB = np.empty(len(smesh))
lABPrev = [0, 0]
for i in range(len(smesh)):
    SSProfile("lAB Rootfinding").tic()
    s = smesh[i]
    lAB = l_a_l_b_rootfinding(s, lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)
    # print(lAB)
    la[i] = lAB[0]
    lb[i] = lAB[1]
    hLALB[i] = lb[i]+la[i]
    lABPrev = lAB
    SSProfile("lAB Rootfinding").toc()
end = time.time()
print("end lAB method")
labtime = end-start
print("lAB time:",end-start)

# Solve the governing equations to get the implied input function values
# This method is a bit incomplete, as it assumes rn agrees where used as an input to the method
# However, since it provides low error, it is probably fine
checkcI = self.t*rn*(lb+la)*(lb/2-la/2)
checkrn = (la+lb)/(np.log((rn+lb)/(rn-la)))
# Check that the input and resultant rn, cI agree
plt.figure("Checking rootfinding accuracy")
plt.plot(smesh, checkcI-PPoly_Eval(smesh, geometryDef[2]), label="cI error")
plt.plot(smesh, checkrn-r_n(smesh, geometryDef[0], geometryDef[1]), label="rn error")
plt.legend()


# Numerically estimate the derivatives of la, lb, rn
numericalDerivatives = np.empty([3, len(smesh)])

for i in range(len(smesh)):
    if i==0:
        numericalDerivatives[0,i]=(la[i+1]-la[i])/(smesh[i+1])
        numericalDerivatives[1,i]=(lb[i+1]-lb[i])/(smesh[i+1])
        numericalDerivatives[2,i]=(rn[i+1]-rn[i])/(smesh[i+1])
    if i==len(smesh)-1:
        numericalDerivatives[0,i]=(la[i]-la[i-1])/(smesh[i]-smesh[i-1])
        numericalDerivatives[1,i]=(lb[i]-lb[i-1])/(smesh[i]-smesh[i-1])
        numericalDerivatives[2,i]=(rn[i]-rn[i-1])/(smesh[i]-smesh[i-1])
    else:
        numericalDerivatives[0,i]=(la[i+1]-la[i-1])/(smesh[i+1]-smesh[i-1])
        numericalDerivatives[1,i]=(lb[i+1]-lb[i-1])/(smesh[i+1]-smesh[i-1])
        numericalDerivatives[2,i]=(rn[i+1]-rn[i-1])/(smesh[i+1]-smesh[i-1])

dlads = numericalDerivatives[0,:]
dlbds = numericalDerivatives[1,:]
drnds = numericalDerivatives[2,:]

# get vectors for the analytical derivatives of different quantities (drnds is Right)
dcIds  = d_cI_d_s(smesh, geometryDef[2])
dads   = d_alpha_d_s(smesh, geometryDef[0], geometryDef[1])
d2ads2 = d_2_alpha_d_s_2(smesh, geometryDef[0], geometryDef[1])
drndsAnalytical  = d_rn_d_s(smesh, geometryDef[0], geometryDef[1])
cI     = PPoly_Eval(smesh, geometryDef[2])

# This proves drnds is wrong
plt.figure("analytical vs numerical rn differentiation")
plt.plot(drnds, label="numerical")
plt.plot(drndsAnalytical, label="analytical")
plt.plot(rn, label = "rn")
plt.legend()

# this also proves its wrong by plotting the error between them
checkdrnds = ((dlads+dlbds)/np.log((rn+lb)/(rn-la))
              -(la+lb)*(1/(np.log((rn+lb)/(rn-la))**2))
              *((rn-la)/(rn+lb))
              *((drnds+dlbds)/(rn-la)-(rn+lb)*(drnds-dlads)/(rn-la)**2))
checkdcIds = (self.t*(drnds*(lb**2/2-la**2/2)+(lb*dlbds-la*dlads)*rn))

plt.figure("rn, cI derivative check")
plt.plot(drnds-checkdrnds, label="drnds error")
plt.plot(dcIds-checkdcIds, label="dcIds error")
plt.legend()

# define geometric functions used in differential equation
geoFunc1 = (dcIds*1/rn-cI*drnds/rn**2)/self.t
geoFunc2 = rn*(-drnds*((la+lb)/rn**2+(-la-lb)/((rn+lb)*(rn-la))))
# # print(geoFunc1, geoFunc2)
geoFuncs = np.array([geoFunc1, geoFunc2])

if False:
    # Use the ODE derived from governing equations to calculate algebraic derivatives for la, lb
    altAlgebraicDerivatives = np.empty([2, len(smesh)])
    for i in range(len(smesh)):
        altAlgebraicDerivatives[0,i], altAlgebraicDerivatives[1,i] = geo_ODE(smesh[i], [la[i], lb[i]], [[geometryDef]])

    plt.figure("compare ODE to numerical derivative")
    plt.plot(smesh, dlads, label="numeric dlads")
    plt.plot(smesh, dlbds, label="numeric dlbds")
    plt.plot(smesh, altAlgebraicDerivatives[0,:], label="algebraic dlads")
    plt.plot(smesh, altAlgebraicDerivatives[1,:], label="algebraic dlbds")
    plt.plot(smesh, dlads-altAlgebraicDerivatives[0,:], label="dlads error")
    plt.plot(smesh, dlbds-altAlgebraicDerivatives[1,:], label="dlbds error")
    plt.legend()

    # Check whether the numerical derivatives of la, lb satisfy the diffEQ
    plt.figure("Do numerical derivatives satisfy diffEQ?")
    plt.plot((geoFunc1-(-la*numericalDerivatives[0,:]+lb*numericalDerivatives[1,:])), label="numeric function 1 (Ic)")
    plt.plot((geoFunc2-((rn/(rn-la)-1)*numericalDerivatives[0,:] + (rn/(rn+lb)-1)*numericalDerivatives[1,:])), label="numeric function 2 (Rn)")
    plt.legend()

    plt.figure("Do numerical derivatives satisfy diffEQ? (Direct Error)")
    plt.plot(smesh, numericalDerivatives[0,:]-altAlgebraicDerivatives[0,:], label="dlads error")
    plt.plot(smesh, numericalDerivatives[1,:]-altAlgebraicDerivatives[1,:], label="dlbds error")
    plt.plot(smesh, np.transpose(numericalDerivatives), label="numerical rootfinding derivatives")
    plt.plot(smesh, np.transpose(altAlgebraicDerivatives), label="ODE algebraic derivatives")
    plt.axhline(0)
    plt.ylim(-0.3, 0.5)
    plt.legend()

# a bunch of bullshit to plot the shape of the spring
ecc = Ic/(self.t*hLALB*rn)
xb = -lb*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
xa = -la*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yb = lb*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
ya = la*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

xrc = ecc*np.sin(alpha_xy(smesh, geometryDef[0], geometryDef[1]))
yrc = ecc*np.cos(alpha_xy(smesh, geometryDef[0], geometryDef[1]))

# print(xrc, yrc)
xorg = coord(smesh, geometryDef[0])
yorg = coord(smesh, geometryDef[1])

plt.figure("geometry results")
plt.plot(xorg+xb,yorg+yb)
plt.plot(xorg-xa,yorg-ya)
plt.plot(-(xorg+xb),-(yorg+yb))
plt.plot(-(xorg-xa),-(yorg-ya))
plt.plot(xorg+xrc, yorg+yrc, label="AB rootfinding centroidal axis")
plt.plot(xorg, yorg, label="rn")
plt.axis('equal')
plt.legend()

# figure out if the numerical derivatives create the same geometry functions
backwardsGeoFunctions = np.empty([2, len(smesh)])
independentCheck = np.empty(len(smesh))
for i in range(len(smesh)):
    states   = np.array([[-la[i], lb[i]], [rn[i]/(rn[i]-la[i])-1, rn[i]/(rn[i]+lb[i])-1]])
    # states   = np.array([[-la[i], lb[i]], [-1/la[i], 1/lb[i]]])
    qdot     = np.array([[numericalDerivatives[0,i]], [numericalDerivatives[0,i]]])
    
    backwardsGeoFunctions[0,i], backwardsGeoFunctions[1,i] = (states.dot(qdot))
    independentCheck[i] = backwardsGeoFunctions[0,i]/backwardsGeoFunctions[1,i]

candidate = (geoFunc1/geoFunc2)

plt.figure("geofunction compare")
plt.plot(smesh, np.transpose(backwardsGeoFunctions), label="apparent functions")
plt.plot(smesh, np.transpose(geoFuncs), label="derived functions")

plt.ylim(-.1, .2)
plt.legend()

print("start forward integration process")
start = time.time()
offset = 1
lABprev = [la[offset], lb[offset]]
lAB0 = l_a_l_b_rootfinding(smesh[offset], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)
print(lAB0)
res    =  fixed_rk4(geo_ODE, lAB0, smesh[offset:-1], geometryDef)
print([[[geometryDef]]])
# print("about to variable mesh solve")
# altRes = variable_mesh_solve(geo_ODE, [smesh[offset], smesh[-1]], lAB0, args=[[[geometryDef]]], method="DOP853", t_eval=smesh[offset:-1])
# print("done variable mesh solving")

laForward = res[0,:]
lbForward = res[1,:]
hLALBF = laForward+lbForward

# laForwardV    = altRes.y[0,:]
# lbForwardV    = altRes.y[1,:]
# variableSmesh = altRes.t
# hLALBFV = laForwardV+lbForwardV

end = time.time()
print("end forward integration method")
forwardTime = end-start
print("forward time", forwardTime)


plt.figure("compare ODE to numerical derivative")
plt.plot(smesh, dlads, label="numeric dlads")
plt.plot(smesh, dlbds, label="numeric dlbds")
plt.plot(np.array(geo_ODE.all_inputs_ever),np.array(geo_ODE.all_outputs_ever)[:,0], label="intro dlads")
plt.plot(np.array(geo_ODE.all_inputs_ever),np.array(geo_ODE.all_outputs_ever)[:,1], label="intro dlbds")
plt.legend()

plt.figure("compare ODE's p inputs to correct inputs")
plt.plot(smesh, la, label="true la")
plt.plot(smesh, lb, label="true lb")
plt.plot(np.array(geo_ODE.all_inputs_ever),np.array(geo_ODE.all_p_inputs_ever)[:,0], label="intro la")
plt.plot(np.array(geo_ODE.all_inputs_ever),np.array(geo_ODE.all_p_inputs_ever)[:,1], label="intro lb")
plt.legend()


plt.figure("direct result compare")
plt.ylim(0, 1)
plt.plot(smesh, la, label='rootfinding la')
plt.plot(smesh, lb, label='rootfinding lb')
plt.plot(smesh[offset:-1], laForward, label='laF fixed mesh')
plt.plot(smesh[offset:-1], lbForward, label='lbF fixed mesh')
# plt.plot(variableSmesh, laForwardV, label='laF variable mesh')
# plt.plot(variableSmesh, lbForwardV, label='lbF variable mesh')
plt.legend()


kn = dads
dknds = d2ads2

states = np.array([la*kn*(1+lb*kn), -lb*kn*(1-la*kn)])
RHS = dknds*(la+lb)*(lb-la-la*lb*kn)
derivatives = np.array([dlads,dlbds])
plt.figure("checking refactor attempt")
plt.plot(states[0,:]*derivatives[0,:]+states[1,:]*derivatives[1,:])
plt.plot( RHS)
plt.plot(dknds)
plt.plot(kn)

plt.show()