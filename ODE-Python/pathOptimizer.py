from os import path
from xml.etree.ElementInclude import default_loader
import numpy as np
import scipy.optimize
from modules.PATHDEF import *
from scipy.optimize import minimize, NonlinearConstraint, differential_evolution
from matplotlib import pyplot as plt

prevArcLen = 6


def end_constraints(designVector, pathObject: RadiallyEndedPolynomial, granuality, originalEndPoint):
    designVector[-1] = 0
    designVector[int(len(designVector)/2-1)] = 1.1
    update_path(designVector, pathObject)
    lastPoint = pathObject.get_neutralSurface(granuality)[-1]
    return np.linalg.norm(lastPoint-originalEndPoint, 2)

def limit_constraints(designVector, pathObject: RadiallyEndedPolynomial):
    designVector[-1] = 0
    designVector[int(len(designVector)/2-1)] = 1.1
    update_path(designVector, pathObject)
    outerLimit = pathObject.outerRadius
    innerLimit = pathObject.innerRadius
    path = pathObject.get_neutralSurface(100)
    totalViolatingPoints = 0
    for i in path:
        if np.linalg.norm(i)<=innerLimit or np.linalg.norm(i)>=outerLimit:
            totalViolatingPoints+=1
    return totalViolatingPoints

def objective(designVector, pathObject: RadiallyEndedPolynomial, granuality=100):
    # print(designVector)
    # print(pathObject)
    designVector[-1] = 0
    designVector[int(len(designVector)/2-1)] = 1.1
    update_path(designVector, pathObject)
    length = pathObject.measure_length()
    mesh = np.linspace(0,length,granuality+1)
    curvature = pathObject.get_rn(mesh)
    curvature = np.linalg.norm(curvature, np.inf)
    return - length

def update_path(designVector, pathObject: RadiallyEndedPolynomial):
    XCoeffs = designVector[0:int(len(designVector)/2)]
    YCoeffs = designVector[int(len(designVector)/2):int(len(designVector))]
    pathObject.XCoeffs = XCoeffs
    pathObject.YCoeffs = YCoeffs

def optimizePath(pathObject: RadiallyEndedPolynomial, designVector):
    constraints = [{"type": "ineq", "fun": end_constraints, "args": (pathObject, 100, pathObject.pts[-1])},
                   {"type": "ineq", "fun": limit_constraints, "args": [pathObject]}]
    # constraints = [NonlinearConstraint(end_constraints, -1e-6, 1e-6), NonlinearConstraint(limit_constraints, 0, 0)]
    res = minimize(objective, designVector, args=pathObject, method='SLSQP')
    return res

# def globalOptimizePath(pathObject: RadiallyEndedPolynomial, designVector):
#     res = differential_evolution()

def main():
    path   = RadiallyEndedPolynomial(1, 6)
    designVector = np.vstack((path.XCoeffs, path.YCoeffs)).flatten()
    defaultCurve = path.get_neutralSurface(100)
    plt.plot(defaultCurve[:,0], defaultCurve[:,1])
    # print(designVector)
    res = optimizePath(path, designVector)
    print(designVector)
    print(res.x)
    print(res.fun)    
    print(res.success)
    print(res.message)
    # print(res.nit)
    optimizedCurve = path.get_neutralSurface(100)
    plt.plot(optimizedCurve[:,0], optimizedCurve[:,1])
    
    plt.axis('equal')
    # plt.axes('equal')
    try:
        plt.show()
    except KeyboardInterrupt:
        return 0

if __name__ == "__main__":
    main()
