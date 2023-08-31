import math
from numpy import number
import numpy as np
from FrequencyCalculations import Frequency_Tides, Frequency_Tides_Masses
from const import *
import scipy as sc
from scipy import integrate, linalg
import matplotlib.pyplot as plt
from HelperCalculations import getParamsWithStep, derivative

def GWsignal(t, t_obs, params):
    #generates simple sinusoidal signal
    A = params[0]
    phi0 = params[1]
    alpha = params[2]
    beta = params[3]
    gamma = params[4]
    #parameterize the signal - now params around order unity
    phi_parameterized = phi0 + 2*np.pi*alpha*(t/t_obs) + (np.pi)*beta*(t/t_obs)**2 + (np.pi/3)*gamma*(t/t_obs)**3
    #phi = phi0 + 2*np.pi*f0*(t) + np.pi*fD*(t)**2 + (np.pi/3)*fDD*(t)**3
    h = A * np.cos(phi_parameterized)
    return h

def GWsignal_masses(t, t_obs, params):
    #generates simple sinusoidal signal
    A = params[0]
    phi0 = params[1]
    f0 = params[2]
    mass1 = params[3]
    mass2 = params[4]

    fD, fDD = Frequency_Tides_Masses(f0, mass1, mass2, t_obs)
    return GWSignal_Tides(t, A, phi0, f0, fD, fDD)

def GWsignal_chirpTotal(t, t_obs, params):
    #generates simple sinusoidal signal
    A = params[0]
    phi0 = params[1]
    f0 = params[2]
    chirp = params[3]
    total = params[4]

    fD, fDD = Frequency_Tides(f0, chirp, total, t_obs)
    return GWSignal_Tides(t, A, phi0, f0, fD, fDD)


def GWSignal_Tides(t, A, phi0, f0, fD, fDD):
    #parameterize the signal - now params around order unity
    phi = phi0 + 2*np.pi*f0*(t) + np.pi*fD*(t)**2 + (np.pi/3)*fDD*(t)**3
    h = A * np.cos(phi)
    return h


def getSNR(t_obs, A, phi0, f0, fD, fDD):
    params = np.array([A, phi0, f0, fD, fDD])
    integrand = lambda t: GWsignal(t, t_obs, params) * GWsignal(t, t_obs, params)
    SNR = integrate.quad(integrand, 0, t_obs)
    return(np.sqrt(SNR[0] * 2))

def getFisherMatrix(t_obs, func, params, labels_params = []):
    # generate labels / assert the number is right
    if len(labels_params) == 0:
        labels_params = np.arange(0, len(params))
    else:
        assert(len(params) == len(labels_params))

    # generate step size
    h_params = []
    for param in params:
        order_of_magnitude = math.floor(math.log10(param))
        order_of_magnitude -= 6
        h_params.append(10**order_of_magnitude)

    #initialize matrix and step size vectors
    fisher = np.zeros((np.size(params), np.size(params)))

    #set up finite difference method - upper and lower step sizes more robust
    get_fx = lambda t : lambda p : func(t, p)

    for i in range(len(params)):
        for j in range(len(params)):
            di = lambda t : derivative(get_fx(t), params, h_params, i)
            dj = lambda t : derivative(get_fx(t), params, h_params, j)
            integrationFunction = lambda t : di(t) * dj(t)
            value = integrate.quad(integrationFunction, 0, t_obs, limit=200)
            fisher[i, j] = value[0] * 2
            
    #print(fisher)
    for i in range(len(params)):
        print('Parameter %s: %e'%(labels_params[i], params[i]))

    sigma = linalg.inv(fisher)
    #print(sigma)

    for i in range(len(params)):
        print('Error in %s: %e'%(labels_params[i], np.sqrt(sigma[i, i])))

    return fisher