import numpy as np
from scipy.integrate import solve_ivp as ODE45
import matplotlib.pyplot as plt

deg2rad = np.pi/180

Rin = 1
Rout = 6
dBMax = 5*deg2rad
thmin = 30*deg2rad
thmax = 360*deg2rad
q = 1.4
c = (Rout-Rin)/thmax**q
h = .25
t = .5
E = 100000


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

def s(theta):
    s = integral(sintegrand, theta)
    return s

def deriv_estimate(fun, theta):
    derivEstimate = (fun(theta+.001)-fun(theta-.001))/.002
    return derivEstimate

def sintegrand(theta):
    sintegrand = np.sqrt(deriv_estimate(xyn, theta)[0]**2+deriv_estimate(xyn, theta)[1]**2)
    return sintegrand

def integral(fun, theta):
    #centered Riemman sum
    res = 1000
    step = (theta-thmin)/res
    i = 1
    sum = 0
    while i < res:
        sum += fun(thmin+i*step)*2*step
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

def deform_ODE(theta, gamma):
    Fx = F*np.sin(alpha(thmax))
    Fy = F*np.cos(alpha(thmax))
    dgSds = deriv_estimate(gen_stiffness, theta)
    LHS = np.empty(2)
    LHS[0] = gamma[1]
    LHS[1] = (Fx*np.sin(alpha(theta)+gamma[0])-Fy*np.cos(alpha(theta)+gamma[0])-dgSds*gamma[1])/gen_stiffness(theta)
    return LHS

print(c)
print(np.linspace(thmin, thmax, 81))
S = np.empty(81)
for i in range(len(np.linspace(thmin, thmax, 81))):
    S[i] = s(np.linspace(thmin, thmax, 81)[i])

# plt.plot(np.linspace(thmin, thmax, 81), S)
print(S[-1])

# F = 50
# th_span = [thmin, thmax]
# smax = s(thmax)
# s_span  = [0, smax]
# gamma = np.zeros(2)
# dfmation = ODE45(deform_ODE, s_span, gamma, dense_output=True, options={'max_step': 1.0})

# plt.plot(dfmation.t, dfmation.sol(dfmation.t)[0])
# print(dfmation.message)
plt.show()
