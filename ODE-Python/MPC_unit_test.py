import time
import math
import numpy as np
import scipy as sp
from scipy import integrate
from scipy import optimize as op
import numpy.linalg as lin
import scipy.linalg as slin
import matplotlib.pyplot as plt
import matplotlib.transforms as tfm
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *
from openpyxl import *

# example code from internet

mu = 0.02
l = 0.0279
eta = 0.01

def fun_measles(x, y):
    beta = 1575 * (1 + np.cos(2 * np.pi * x))
    return np.vstack((
        mu - beta * y[0] * y[2],
        beta * y[0] * y[2] - y[1] / l,
        y[1] / l - y[2] / eta,

    ))

def bc_measles(ya, yb):
    return ya - yb

x_measles = np.linspace(0, 1, 100)
y_measles = np.full((3, x_measles.shape[0]), 0.01)

res_measles = integrate.solve_bvp(fun_measles, bc_measles, x_measles, y_measles, verbose=2)

x_measles_plot = np.linspace(0, 1, 100)
y_measles_plot = res_measles.sol(x_measles_plot)