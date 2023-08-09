from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg
import matplotlib.pyplot as plt
import collections.abc

def getChirpMass(mass1, mass2):
    return((mass1*mass2)**(3/5) / (mass1 + mass2)**(1/5))

def getTotalMass(mass1, mass2):
    return(mass1 + mass2)

def derivative(fx, params, h, i = 0):
    isParamsArray = isinstance(params, collections.abc.Sequence)
    isHArray = isinstance(h, collections.abc.Sequence)

    if isParamsArray != isHArray:
        print("ERROR: Params and h must be the same size")
        return "ERROR"
    elif isParamsArray != True:
        params = [params]
        h = [h]
    elif len(params) != len(h):
        print("ERROR: Params and h must be the same size")
        return "ERROR"

    return (fx(getParamsWithStep(params, i, h[i])) - fx(getParamsWithStep(params, i, h[i], False))) / (2*h[i])

def getParamsWithStep(params, target, step, stepUp = True):
    newParams = np.copy(params)

    if (stepUp):
        newParams[target] = params[target] + step
    else:
        newParams[target] = params[target] - step

    return newParams