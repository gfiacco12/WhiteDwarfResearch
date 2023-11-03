import numpy as np
from const import *

def getChirpMass(mass1, mass2):
    return(((mass1*mass2)**(3/5)) / ((mass1 + mass2)**(1/5)))

def getTotalMass(mass1, mass2):
    return(mass1 + mass2)

def convertChirpTotal_to_M1M2(chirpmass, totalmass):
    m1 = (1/2) * (totalmass + np.sqrt(totalmass**2 - (4*chirpmass**(5/3)*totalmass**(1/3))))
    m2 = (1/2) * (totalmass - np.sqrt(totalmass**(1/3) * (totalmass**(5/3) - 4*chirpmass**(5/3))))
    return [m1,m2]

def get_comp_mass_Mc_Mt(Mc, Mt):
    eta = (Mc/Mt)**(5./3.)
    qq = ((1-2.*eta) - np.sqrt((1-2.*eta)**2. - 4.*eta**2.))/(2.*eta)
    M1=Mt/(1.+qq)
    M2=qq*Mt/(1.+qq)
    return M1, M2 

def calculateAmplitude(SNR, t_obs):
    #calculates approx amplitude from SNR^2 = int( h(t)^2 dt) - time avg over 0 to t_obs
    return(SNR / np.sqrt(t_obs))

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