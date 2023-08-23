from numpy import number
import numpy as np
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

def getSNR(t_obs, A, phi0, f0, fD, fDD):
    params = np.array([A, phi0, f0, fD, fDD])
    integrand = lambda t: GWsignal(t, t_obs, params) * GWsignal(t, t_obs, params)
    SNR = integrate.quad(integrand, 0, t_obs)
    return(np.sqrt(SNR[0] * 2))

def getFisherMatrix(t_obs, A, phi0, f0, fD, fDD):
    #define parameters
    alpha = f0 * t_obs
    beta = fD * (t_obs**2)
    gamma = fDD * (t_obs**3)
    params = np.array([A, phi0, alpha, beta, gamma])
    h_param = np.array([1.e-6, 0.001, 1.e1, 0.1, 0.0001])
    
    label = [r'$A$', r'$\phi_{0}$',  r'$\alpha}$', r'$\beta$', r'$\gamma$']
    
    """ h = np.array([1.e-6, 0.001, 1.e-6, 1.e-17, 1.e-29])
    label2 = [r'$A$', r'$\phi_{0}$',  r'$f_{0}$', r'$\dot{f}$', r'$\ddot{f}$'] """

    #initialize matrix and step size vectors
    fisher = np.zeros((np.size(params), np.size(params)))

    #set up finite difference method - upper and lower step sizes more robust
    get_fx = lambda t : lambda p : (GWsignal(t, t_obs, p))

    for i in range(len(params)):
        for j in range(len(params)):
            di = lambda t : derivative(get_fx(t), params, h_param, i)
            dj = lambda t : derivative(get_fx(t), params, h_param, j)
            integrationFunction = lambda t : di(t) * dj(t)
            value = integrate.quad(integrationFunction, 0, t_obs, limit=200)
            fisher[i, j] = value[0] * 2
    print(fisher)
    for i in range(len(params)):
        print('Parameter %s: %e'%(label[i], params[i]))
    sigma = linalg.inv(fisher)
    print("alp:", f0*t_obs)
    print("beta:", fD * (t_obs**2))
    print("gamma:", fDD * (t_obs**3))

    for i in range(len(params)):
        print('Error in %s: %e'%(label[i], np.sqrt(sigma[i, i])))
    return fisher