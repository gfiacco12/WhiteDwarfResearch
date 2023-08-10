import numpy as np
from const import *

def getChirpMass(mass1, mass2):
    return(((mass1*mass2)**(3/5)) / ((mass1 + mass2)**(1/5)))

def getTotalMass(mass1, mass2):
    return(mass1 + mass2)

def calculateAmplitude(SNR, t_obs):
    #calculates approx amplitude from SNR^2 = int( h(t)^2 dt) - time avg over 0 to t_obs
    return(np.sqrt(2)*SNR / np.sqrt(t_obs))

def calculateAmplitude_phys(dl, chirpMass, freq0):
    #calculates approx amplitude from SNR^2 = int( h(t)^2 dt) - time avg over 0 to t_obs
    return(np.pi**2/3 * chirpMass**(5/3) * freq0**(2/3) / dl)

def derivative(fx, params, h, i = 0):
    np.isscalar(params)
    isParamsArray = not np.isscalar(params)
    isHArray = not np.isscalar(h)

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