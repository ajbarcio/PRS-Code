import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
from materials import *

#def forward_sum(array):
#    ret = np.array(array)
#    for i in range(len(ret))[1:]:
#        ret[i]+=ret[i-1]
#    return ret

#x = np.linspace(0,1,21)
# print(x)
#for i in range(len(x)):
    #y = x

#integ = forward_sum(y)
#print(integ)

A = np.array([3])
print(A.searchsorted(3))