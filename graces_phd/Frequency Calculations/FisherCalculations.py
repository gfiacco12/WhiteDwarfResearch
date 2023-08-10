from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg
import matplotlib.pyplot as plt
from HelperCalculations import getParamsWithStep

def GWsignal(t, t_obs, params):
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

    phi = phi0 + 2*np.pi*alpha*(t/t_obs) + np.pi*beta*(t/t_obs)**2 + (np.pi/3)*gamma*(t/t_obs)**3
    h = A * np.cos(2 * np.pi * phi)
    return h

def getSNR(t_obs, A, phi0, f0, fD, fDD):
    params = np.array([A, phi0, f0, fD, fDD])
    integrand = lambda t: GWsignal(t, t_obs, params) * GWsignal(t, t_obs, params)
    SNR = 2 * integrate.quad(integrand, 0, t_obs)
    return(np.sqrt(SNR[0]))

def getFisherMatrix(t_obs, A, phi0, f0, fD, fDD):
    #define parameters
    params = np.array([A, phi0, f0, fD, fDD])
    h = np.array([[1.e-6, 0.001, 1.e3, 0.1, 0.001],
                    [-1.e-6, -0.001, -1.e3, -0.1, -0.001]])
    label = [r'$A$', r'$\phi_{0}$',  r'$f_{0}$', r'$\dot{f}_{0}$', r'$\ddot{f}_{0}$']

    #initialize matrix and step size vectors
    fisher = np.zeros((np.size(params), np.size(params)))

    #set up finite difference method - upper and lower step sizes more robust
    """ fx = lambda x, t : (GWsignal(t, t_obs, getParamsWithStep(params, x, h[x], True)) - 
                        GWsignal(t, t_obs, getParamsWithStep(params, x, h[x], False))) / (2 * h[x])
    """
    #hang's method:  
    for i in range(len(params)):
        par_u = params.copy()
        par_u += h[0,:]
        par_l = params.copy()
        par_l += h[1,:]

        fx = lambda x, t: (GWsignal(t, t_obs, par_u) -  GWsignal(t, t_obs, par_l)) / (h[0, x] - h[1, x])
        for j in range(len(params)):
            integrationFunction = lambda t : fx(i, t) * fx(j, t)
            value = 2 * integrate.quad(integrationFunction, 0, t_obs)
            fisher[i][j] = value[0]
    #print(par_l)
    #print(par_u)
    #print(fisher)
    #print("---------")
    sigma = linalg.inv(fisher)

    for i in range(len(params)):
        print('Error in %s: %e'%(label[i], np.sqrt(sigma[i, i])))
    return fisher

#From Hangs code
def inner_product(h1, h2, freq, psd):  
    """
    inner product
    """
    integrand = (np.conj(h1)*h2+h1*np.conj(h2))/psd
    rho_sq = 2.*integrate.trapz((integrand/2), freq)
    return rho_sq

def fisher(t_obs, A, phi0, f0, fD, fDD, psd):
    """
    numerically compute the fisher matrix
    """
    t = np.linspace(0, t_obs, 1000)
    params = np.array([A, phi0, f0, fD, fDD])
    nPt=len(t)
    nDof=len(params)
    dpar= np.array([[1.e-6, 0.001, 1.e3, 0.1, 0.001],
                    [-1.e-6, -0.001, -1.e3, -0.1, -0.001]])
    dh=np.zeros([nDof, nPt], dtype=np.complex128)
    for i in range(nDof):
        par_u = params.copy()
        par_u[i] += dpar[0, i]
        hh_u = GWsignal(t, t_obs, par_u)
        
        par_l = params.copy()
        par_l[i] += dpar[1, i]
        hh_l = GWsignal(t, t_obs, par_l)
        
        dh[i, :] = (hh_u - hh_l) / (dpar[0, i] - dpar[1, i])
        
    gamma=np.zeros([nDof,nDof])
    for i in range(nDof):
        for j in range(i, nDof, 1):
            gamma[i, j]=np.real(inner_product(dh[i, :], dh[j, :], t, psd))
    
    for i in range(nDof):
        for j in range(i):
            gamma[i, j]=gamma[j, i]
    print(gamma)
    print("---------")
    sigma = linalg.inv(gamma)
    label = [r'$A$', r'$\phi_{0}$', r'$\alpha$', r'$\beta$', r'$\gamma$']
    for i in range(len(params)):
        print('Error in %s: %e'%(label[i], np.sqrt(sigma[i, i])))
            
    return gamma