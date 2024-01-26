import numpy as np
from scipy.integrate import solve_ivp as ODE45
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

deg2rad = np.pi/180
rad2deg = 1/deg2rad

Rin = 1
Rout = 6
dBMax = 5*deg2rad
thmin = 30*deg2rad
thmax = 360*deg2rad
q = 1.4
c = (Rout-Rin)/thmax**q
h = 0.5
hroot = 0.5
hmin = 0.25
t = .5
E = 100000000


def R(theta):
    R = c*theta**q + Rin
    return R

def R_prime(theta):
    RPrime = (c*q*theta**(q-1))
    return RPrime

def R_prime_prime(theta):
    RPrimePrime = c*q*(q-1)*theta**(q-2)
    return RPrimePrime

def r_c(theta):
    Rc = ((R(theta)**2+R_prime(theta)**2)**(3/2))/ \
         np.abs(R(theta)**2+2*R_prime(theta)**2-R(theta)*R_prime_prime(theta))
    return Rc
    
def r_n(theta):
    Rn = h/(np.log((r_c(theta)+.5*h)/(r_c(theta)-.5*h)))
    return Rn

def e(theta):
    e = r_c(theta)-r_n(theta)
    return e

def xyc(theta):
    xyc = np.array([R(theta)*np.cos(theta), R(theta)*np.sin(theta)])
    return xyc

def xyn(theta):
    xn = xyc(theta)[0] - e(theta)*R(theta)*np.cos(theta)/ \
         np.sqrt((R(theta)*-np.sin(theta))**2+(R(theta)*np.cos(theta))**2)
    yn = xyc(theta)[1] - e(theta)*R(theta)*-np.sin(theta)/ \
         np.sqrt((R(theta)*-np.sin(theta))**2+(R(theta)*np.cos(theta))**2)
    xyn = np.array([xn, yn])
    return xyn

def s_exact(theta):
    s = integral_theta(sintegrand, theta)
    return s

def deriv_estimate(fun, theta):
    derivEstimate = (fun(theta+.001)-fun(theta-.001))/.002
    return derivEstimate

def sintegrand(theta):
    sintegrand = np.sqrt(deriv_estimate(xyn, theta)[0]**2+deriv_estimate(xyn, theta)[1]**2)
    return sintegrand

def integral_theta(fun, theta):
    #centered Riemman sum
    res = 1000
    step = (theta-thmin)/res
    i = 1
    sum = 0
    while i < res:
        sum += fun(thmin+i*step)*2*step
        i += 2
    return sum

def integral_s(fun, s):
    #centered Riemman sum
    res = 1000
    step = (s)/res
    i = 1
    sum = 0
    while i < res:
        sum += fun(i*step)*2*step
        i += 2
    return sum

def alpha(theta):
    alpha = np.arctan2(deriv_estimate(xyn, theta)[1],deriv_estimate(xyn, theta)[0])
    return alpha

def gen_stiffness(theta):
    if np.isinf(r_n(theta)) and e(theta) == 0:
        genStiffness = t*h*E
    else:
        genStiffness = t*h*E*e(theta)*r_n(theta)
    return genStiffness

def stupid_cos(s):
    return np.cos(alpha(theta_of_s_approx(s, *coeffs))+gamma_fun(s))

def stupid_sin(s):
    return np.sin(alpha(theta_of_s_approx(s, *coeffs))+gamma_fun(s))

def xyc_deformed(s):
    xyc = np.empty(2)
    xyc[0] =  integral_s(stupid_cos, s)
    xyc[1] =  integral_s(stupid_sin, s)

def deform_ODE(s, gamma):
    Fx = F*np.sin(alpha(thmax))
    Fy = F*np.cos(alpha(thmax))
    theta = theta_of_s_approx(s, *coeffs)
    dgSds = deriv_estimate(gen_stiffness, theta)/deriv_estimate(s_exact, theta)
    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha(theta)+gamma[0])-Fy*np.cos(alpha(theta)+gamma[0])-dgSds*gamma[1])/gen_stiffness(theta)
    return LHS

#print(c)
# #print(np.linspace(thmin, thmax, 81))
S = np.empty(81)
for i in range(len(np.linspace(thmin, thmax, 81))):
    S[i] = s_exact(np.linspace(thmin, thmax, 81)[i])

#print(S[-1])
# print(np.column_stack([np.linspace(thmin, thmax, 81),S]))

def s_of_theta_approx(x, a, qn):
    return (x-a)**qn-(thmin-a)**qn

def theta_of_s_approx(x, a, qn):
    return a + (x+(thmin-a)**qn)**(1./qn)

coeffs, shit = curve_fit(s_of_theta_approx,np.linspace(thmin, thmax, 81),S,p0=[.52,1.695])
# print(coeffs)

# plt.plot((np.linspace(thmin, thmax, 81)), S)
# plt.plot((np.linspace(thmin, thmax, 81)), s_of_theta_approx(np.linspace(thmin, thmax, 81),*coeffs))
# plt.show()

Fdeforms={}
Fmeshes={}
Fangles={}
Stress={}
maxStress={}
Frange = 5
MaxF = 50
for i in range (-Frange,Frange+1):
    if i==-Frange:
        print("entered loop")
    key=str("dfrm_F_"+str(i))

    F = -i*MaxF/Frange
    th_span = [thmin, thmax]
    smax = s_exact(thmax)
    s_span  = [0, smax]
    gamma = np.zeros(2)
    print("thinking about what happens when F =",F)
    tempDeform = ODE45(deform_ODE, s_span, gamma, dense_output=True, method='LSODA')
    Fdeforms[key] = tempDeform.y[0,:]
    Fmeshes[key]  = tempDeform.t
    Fangles[key] = tempDeform.y[0,-1]
    print(Fangles[key])
    Stress[key] = np.abs(E*(1-r_n(theta_of_s_approx(tempDeform.t, *coeffs))/(r_c(theta_of_s_approx(tempDeform.t, *coeffs))-h/2)) \
                    *r_n(theta_of_s_approx(tempDeform.t, *coeffs))*tempDeform.y[1,:])
    maxStress[key] = np.max(Stress[key])

plt.figure(1)
for key in Fdeforms:
    plt.plot(Fmeshes[key], Fdeforms[key])
i = 0
finalAngles = np.empty(len(range(-Frange,Frange+1)))
for key in Fangles:
    finalAngles[i]=Fangles[key]*rad2deg
    i+=1
i = 0
maxStresses = np.empty(len(range(-Frange,Frange+1)))
for key in Fangles:
    maxStresses[i]=maxStress[key]
    i+=1
plt.figure(2)
plt.plot(np.linspace(-Frange,Frange,len(range(-Frange,Frange+1))),finalAngles)
plt.figure(3)
plt.plot(np.linspace(-Frange,Frange,len(range(-Frange,Frange+1))),(maxStresses))


# for key in Fmeshes:
#     print(key)
# for key in Fdeforms:
#     print(key)

def gamma_fun(s):
    return dfmation.sol(s)[0]

# deformed_shape = np.empty([2,len(dfmation.t)])
# for i in range(len(dfmation.t)):
#     deformed_shape[:,i] = xyc_deformed(dfmation.t[i])

# print(deformed_shape)

plt.show()
