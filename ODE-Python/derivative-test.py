import numpy as np
from scipy.integrate import solve_ivp as ODE45
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def wow(theta):
    return np.array([2*theta, theta-1])


def deriv_estimate(fun, theta):
    derivEstimate = (fun(theta+.001)-fun(theta-.001))/.002
    return derivEstimate

fuck = deriv_estimate(wow, 2)
print(fuck)
