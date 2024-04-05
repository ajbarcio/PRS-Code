import numpy as np
from scipy.integrate import solve_ivp as variable_mesh_solve
import numpy.linalg as lin
import matplotlib.pyplot as plt
import time
from stiffness_library import *
from StatProfiler import SSProfile


fullArcLength = 4.1
deg2rad = np.pi/180
dragVector0 = [1, 2, 2.7, 4, \
               5*deg2rad, 10*deg2rad, 15*deg2rad,  \
               .001, .0005, .001, \
               fullArcLength/2,  fullArcLength/3, 2*fullArcLength/3, fullArcLength]
# dragVector0 = [1.0, 2.0, 6.0, 10.0, 1*deg2rad, 2*deg2rad, 3*deg2rad, 0.5, 0.5, 0.5, \
            #    6, 4, 8, 12]
print("initial initial guess", dragVector0)

geometryDef, smesh = drag_vector_spring(dragVector0)
smesh = np.linspace(0, fullArcLength, 5005)
rn = np.empty(len(smesh))
Ic = np.empty(len(smesh))
for i in range(len(smesh)):
    rn[i] = r_n(smesh[i], geometryDef[0], geometryDef[1])
    Ic[i] = cI_s(smesh[i], geometryDef[2])

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

checkcI = outPlaneThickness*rn*(lb+la)*(lb/2-la/2)
checkrn = (la+lb)/(np.log((rn+lb)/(rn-la)))
plt.figure("Checking rootfinding accuracy")
plt.plot(smesh, checkcI-cI_s(smesh, geometryDef[2]), label="cI error")
plt.plot(smesh, checkrn-r_n(smesh, geometryDef[0], geometryDef[1]), label="rn error")
plt.legend()

numericalDerivatives = np.empty([2, len(smesh)])

for i in range(len(smesh)):
    if i==0:
        numericalDerivatives[0,i]=(la[i+1]-la[i])/(smesh[i+1])
        numericalDerivatives[1,i]=(lb[i+1]-lb[i])/(smesh[i+1])
    if i==len(smesh)-1:
        numericalDerivatives[0,i]=(la[i]-la[i-1])/(smesh[i]-smesh[i-1])
        numericalDerivatives[1,i]=(lb[i]-lb[i-1])/(smesh[i]-smesh[i-1])
    else:
        numericalDerivatives[0,i]=(la[i+1]-la[i-1])/(smesh[i+1]-smesh[i-1])
        numericalDerivatives[1,i]=(lb[i+1]-lb[i-1])/(smesh[i+1]-smesh[i-1])

dcIds  = d_cI_d_s(smesh, geometryDef[2])
dads   = d_alpha_d_s(smesh, geometryDef[0], geometryDef[1])
d2ads2 = d_2_alpha_d_s_2(smesh, geometryDef[0], geometryDef[1])
drnds  = d_rn_d_s(smesh, geometryDef[0], geometryDef[1])
cI     = cI_s(smesh, geometryDef[2])

geoFunc1 = (dcIds*1/rn-cI*drnds/rn**2)/outPlaneThickness
geoFunc2 = rn*(-drnds*((la+lb)/rn**2+(-la-lb)/((rn+lb)*(rn-la))))
# # print(geoFunc1, geoFunc2)
geoFuncs = np.array([geoFunc1, geoFunc2])

altAlgebraicDerivatives = np.empty([2, len(smesh)])
for i in range(len(smesh)):
    altAlgebraicDerivatives[0,i], altAlgebraicDerivatives[1,i] = geo_ODE(smesh[i], [la[i], lb[i]], [[geometryDef]])

plt.figure("check derivatives")
plt.plot(geoFunc1-(-la*altAlgebraicDerivatives[0,:]+lb*altAlgebraicDerivatives[1,:]))
plt.plot(geoFunc2-((rn/(rn-la)-1)*altAlgebraicDerivatives[0,:] + (rn/(rn+lb)-1)*altAlgebraicDerivatives[1,:]))
# THE ROOTFINDING-BASED VALUES SATISFY THE ALGEBRAIC FOMRULATION OF THE DERIVATIVES (THIS IS CIRCULAR BUT IT PROVES OUR ALGEBRA IS INTERNALLY CONSISTEN)
plt.plot((geoFunc1-(-la*numericalDerivatives[0,:]+lb*numericalDerivatives[1,:]))/drnds)
plt.plot((geoFunc2-((rn/(rn-la)-1)*numericalDerivatives[0,:] + (rn/(rn+lb)-1)*numericalDerivatives[1,:]))/drnds)

# BUT THE NUMERICAL DERIVATIVES DO NOT SATISFY THE ALGEBRAIC RELATIONS

# ARE THE NUMERICAL DERIVATIVES SELF CONSISTENT (DO THEY INTEGRATE BACK UP TO THE CORRECT VALUES)?
step = smesh[1]-smesh[0]
laIntegrated = np.empty(len(smesh))
lbIntegrated = np.empty(len(smesh))
laIntegrated[0]  = la[0]
lbIntegrated[0]  = lb[0]
for i in range(len(smesh[1:-1])):
    laIntegrated[i+1] = laIntegrated[i] + numericalDerivatives[0,i+1]*step
    lbIntegrated[i+1] = lbIntegrated[i] + numericalDerivatives[0,i+1]*step
plt.figure("integration error")
plt.plot(la-laIntegrated)
plt.plot(lb-lbIntegrated)

rnIntegrated = np.empty(len(smesh))
cIIntegrated = np.empty(len(smesh))
rnIntegrated[0]  = r_n(0, geometryDef[0], geometryDef[1])
cIIntegrated[0]  = cI_s(0, geometryDef[2])
for i in range(len(smesh[1:-1])):
    rnIntegrated[i+1] = rnIntegrated[i] + d_rn_d_s(smesh[i+1], geometryDef[0], geometryDef[1])*step
    lbIntegrated[i+1] = lbIntegrated[i] + d_cI_d_s(smesh[i+1], geometryDef[2])*step
plt.figure("analytical integration error")
plt.plot(rn-rnIntegrated)
plt.plot(cI-cIIntegrated)
plt.figure("rn, cI")
plt.plot(smesh, rn)
plt.plot(smesh, cI)


ecc = Ic/(outPlaneThickness*hLALB*rn)
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

plt.show()
assert(False)

backwardsGeoFunctions = np.empty([2, len(smesh)])
independentCheck = np.empty(len(smesh))
for i in range(len(smesh)):
    states   = np.array([[-la[i], lb[i]], [rn[i]/(rn[i]-la[i])-1, rn[i]/(rn[i]+lb[i])-1]])
    # states   = np.array([[-la[i], lb[i]], [-1/la[i], 1/lb[i]]])
    qdot     = np.array([[numericalDerivatives[0,i]], [numericalDerivatives[0,i]]])
    
    backwardsGeoFunctions[0,i], backwardsGeoFunctions[1,i] = (states.dot(qdot))
    independentCheck[i] = backwardsGeoFunctions[0,i]/backwardsGeoFunctions[1,i]

candidate = (geoFunc1/geoFunc2)

plt.figure(1000)
plt.plot(smesh, np.transpose(backwardsGeoFunctions), label="apparent functions")
plt.plot(smesh, np.transpose(geoFuncs), label="derived functions")
# plt.plot(smesh, la, label = "la")
# plt.plot(smesh, lb, label = "lb")
plt.ylim(-.1, .2)
plt.legend()
# plt.figure(1001)
# plt.plot(smesh, independentCheck, label="ratio between backwards geofunctions")
# # plt.plot(smesh, -hLALB/independentCheck, label = "ratio between thickness and geofunctions")
# # plt.plot(smesh, -hLALB/6, label="negative thickness")
# plt.plot(smesh, rn, label="rn/smesh")
# plt.ylim(-3, 3)
# plt.legend()
# plt.plot(smesh, candidate)

plt.figure("derivative compare")
plt.plot(smesh, np.transpose(numericalDerivatives), label="numerical rootfinding derivatives")
# plt.plot(smesh, np.transpose(algebraicDerivatives), label="calculated algebraic derivatives")
plt.plot(smesh, np.transpose(altAlgebraicDerivatives), label="ODE algebraic derivatives")
plt.axhline(0)
plt.ylim(-0.3, 0.5)
plt.legend()


print("start forward integration process")
start = time.time()
lABprev = [0, 0]
lAB0 = l_a_l_b_rootfinding(smesh[1], lABPrev, geometryDef[0], geometryDef[1], geometryDef[2], False)
# lAB0 = np.array([lAB0, lAB0])

res    =           fixed_rk4(geo_ODE, lAB0, smesh[1:-1], geometryDef)
print([[[geometryDef]]])
# print("about to variable mesh solve")

# altRes = variable_mesh_solve(geo_ODE, [smesh[1], smesh[-1]], lAB0, args=[[[geometryDef]]])

# print("done variable mesh solving")
# print(res)
laForward = res[0,:]
lbForward = res[1,:]
hLALBF = laForward+lbForward
end = time.time()
print("end forward integration method")
forwardTime = end-start
print("forward time", forwardTime)
print(altAlgebraicDerivatives[0,0], altAlgebraicDerivatives[1,0])
print(altAlgebraicDerivatives[0,1]-altAlgebraicDerivatives[1,1], numericalDerivatives[0,1]-numericalDerivatives[1,1])


# plt.figure(0)
# plt.ylim(0, 1)
# plt.plot(smesh, hAB, label='old way')
# plt.plot(smesh, hLALB, label='new way')
# plt.plot(smesh[0:-1], hLALBF, label='forward integration way')
# plt.legend()
plt.figure("direct result compare")
plt.ylim(0, 1)
plt.plot(smesh, la, label='rootfinding la')
plt.plot(smesh, lb, label='rootfinding lb')
plt.plot(smesh[1:-1], laForward, label='laF')
plt.plot(smesh[1:-1], lbForward, label='lbF')
# plt.plot(altRes.t, altRes.y, label="variable size attempt")
# plt.legend()
# plt.figure(1)
# plt.ylim(0, 1)
# plt.plot(smesh, (a+b)/2, label="centroidal radius, ab method")
# plt.plot(smesh, (2*rn-la+lb)/2, label="centroidal radius, lab method")
# plt.plot(smesh[0:-1], (2*rn[0:-1]-laForward+lbForward)/2, label="centroidal radius, forward method")
# plt.plot(smesh, rn, label="neutral radius")
# plt.legend()
plt.show()