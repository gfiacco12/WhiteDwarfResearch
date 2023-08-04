from numpy import number
import numpy as np
from const import *

def Frequency_1PN(freq0, mass1, mass2, dl, t_obs):

    chirpMass = (mass1*mass2)**(3/5) / (mass1 + mass2)**(1/5)
    totalMass = mass1 + mass2
    eta = (chirpMass/totalMass)**(5/3)
    #0PN point particle freqD
    fdot_pp = 96/5*np.pi**(8/3)*freq0**(11/3)*chirpMass**(5/3)
    fd_corr = ((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3)
    #1PN freqD and freqDD
    fdot = fdot_pp * (1 + ((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fddot = fdot_pp * (fdot/freq0) * ((11/3) + (13/3)*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fddot_0PN = (11/3)*(fdot ** 2) / freq0
    delta_fddot = fddot - fddot_0PN 

    amp = np.pi**2/3 * chirpMass**(5/3) * freq0**(2/3) / dl

    #frequency bin
    df = (1/t_obs)
    #df_1 = 3.1e-8 * (1/t_obs)

    print("1PN CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("1PN Freq Derivative =", fdot, "Hz")
    print("1PN Freq Second Derivative =", fddot, "Hz")
    print("1PN fdot corr:", fd_corr)
    print("The 1PN Fdd correction is:", delta_fddot, "Hz")
    print("Frequency bin:", df, "Hz")
    print("amplitude:", amp)
    print("-----------------------------------------------")
    return fdot, fddot


def Frequency_Tides(freq0, mass1, mass2, dl, t_obs):

    chirpMass = (mass1*mass2)**(3/5) / (mass1 + mass2)**(1/5)
    totalMass = mass1 + mass2
    eta = (chirpMass/totalMass)**(5/3)
    #0PN point particle freqD
    fdot_pp = 96/5*np.pi**(8/3)*freq0**(11/3)*chirpMass**(5/3)
    
    I_wd = 8.51e-10 * ((mass1/0.6*MSOLAR)**(1/3) + (mass2/0.6*MSOLAR)**(1/3))
    I_orb = chirpMass**(5/3) / ((np.pi*freq0)**(4/3))
    #1PN freqD and freqDD
    fdot = fdot_pp * (1 + ((3*I_wd/I_orb)/(1 - (3*I_wd/I_orb))) )
    fddot = (11/3)*(fdot_pp * fdot/freq0 )* ((1 - (21/11)*(I_wd/I_orb)) / ((1 - (3*I_wd/I_orb))**2))
    fddot_0PN = (11/3)*(fdot ** 2) / freq0
    delta_fddot = fddot - fddot_0PN 
    fd_corr = (3*I_wd/I_orb)/(1 - (3*I_wd/I_orb))

    amp = np.pi**2/3 * chirpMass**(5/3) * freq0**(2/3) / dl


    #frequency bin
    df = (1/t_obs)
    #df_1 = 3.1e-8 * (1/t_obs)

    print("TIDAL CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("Tides Freq Derivative =", fdot, "Hz")
    print("Tides Freq Second Derivative =", fddot, "Hz")
    print("Tidal fdot corr:", fd_corr)
    print("The Tidal Fdd correction is:", delta_fddot, "Hz")
    print("Frequency bin:", df, "Hz")
    print("-----------------------------------------------")
    return fdot, fddot


def GWsignal(t_obs, params):
    #generates simple sinusoidal signal
    A = params[0]
    phi0 = params[1]
    f0 = params[2]
    fD = params[3]
    fDD = params[4]
    #parameterize the signal - now params around order unity
    alpha = f0 * t_obs
    beta = fD * t_obs**2
    gamma = fDD * t_obs**3
    t = np.linspace(0, t_obs, 1000)

    """ print("alpha:", alpha) #~ 1.2e6
    print("beta:", beta) # ~210
    print('gamma:', gamma) # ~ 0.1 """
    phi = phi0 + 2*np.pi*alpha*(t/t_obs) + np.pi*beta*(t/t_obs)**2 + (np.pi/3)*gamma*(t/t_obs)**3
    h = A * np.cos(2 * np.pi * phi)
    return h


def calculateAmplitude(SNR, t_obs):
    #calculates approx amplitude from SNR^2 = int( h(t)^2 dt) - time avg over 0 to t_obs
    return(np.sqrt(2)*SNR / np.sqrt(t_obs))


def getFisherMatrix(t_obs, A, phi0, f0, fD, fDD):
    #define parameters
    params = np.array([A, phi0, f0, fD, fDD])
    h = np.array([1.e-6, 0.001, 1.e3, 0.1, 0.001])
    
    #initialize matrix and step size vectors
    fisher = np.zeros((np.size(params), np.size(params)))

    fx = lambda x : (GWsignal(t_obs, getParamsWithStep(params, x, h[x], True))[-1] - GWsignal(t_obs, getParamsWithStep(params, x, h[x], False)))[-1] / (2*h[x])

    for i in range(len(params)):
        for j in range(len(params)):
            fisher[i][j] = fx(i) * fx(j)

    return fisher

def getParamsWithStep(params, target, step, stepUp = True):
    newParams = np.copy(params)

    if (stepUp):
        newParams[target] = params[target] + step
    else:
        newParams[target] = params[target] - step

    return newParams