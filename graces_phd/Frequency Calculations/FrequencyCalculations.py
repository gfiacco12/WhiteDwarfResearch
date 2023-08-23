from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg, optimize
import matplotlib.pyplot as plt
from HelperCalculations import *
from FisherCalculations import getParamsWithStep
import sympy as sp
import torch

def Frequency_1PN(freq0, mass1, mass2, dl, t_obs):

    chirpMass = getChirpMass(mass1, mass2)
    totalMass = getTotalMass(mass1, mass2)
    eta = (chirpMass/totalMass)**(5/3)
    #0PN point particle freqD
    fdot_pp = 96/5*np.pi**(8/3)*freq0**(11/3)*chirpMass**(5/3)
    fd_corr = fdot_pp*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3)
    #1PN freqD and freqDD
    ### Fix this ####
    fdot = fdot_pp * (1 + ((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fddot = fdot_pp * (fdot/freq0) * ((11/3) + (13/3)*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fddot_0PN = (11/3)*(fdot ** 2) / freq0
    fddot_1PN = fddot_0PN*( (3/11)*(13/3)*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    delta_fddot = fddot - fddot_0PN 

    amp = np.pi**2/3 * chirpMass**(5/3) * freq0**(2/3) / dl

    #frequency bin
    df = (1/t_obs)
    dfd = (1/t_obs**2)
    dfdd = (1/ t_obs**3)

    print("1PN CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp * t_obs**2, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("1PN Freq Derivative =", fdot, "Hz")
    print("1PN Freq Second Derivative =", fddot, "Hz")
    print("1PN Freq Derivative Dimensionless=", fdot * t_obs**2, "Hz")
    print("1PN Freq Second Derivative Dimensionless =", fddot * t_obs**3, "Hz")
    print("1PN fdot corr:", fd_corr * t_obs**2)
    print("The 1PN Fdd correction is:", delta_fddot * t_obs**3, "Hz")
    print("Test 1PN fdd corr:", fddot_1PN*t_obs**3)
    print("amplitude:", amp)
    print("-----------------------------------------------")
    print("Frequency bin:", df)
    print("Change in Freq bin due to fdot:", dfd)
    print("Change in Freq bin due to fddot:", dfdd)
    print("-----------------------------------------------")
    return fdot, fddot


def Frequency_Tides(freq0, mass1, mass2, dl, t_obs):

    chirpMass = getChirpMass(mass1, mass2)
    totalMass = getTotalMass(mass1, mass2)
    eta = (chirpMass/totalMass)**(5/3)
    #0PN point particle freqD
    fdot_pp = 96/5*np.pi**(8/3)*freq0**(11/3)*chirpMass**(5/3)
    
    I_wd = 8.51e-10 * ((mass1/(0.6*MSOLAR))**(1/3) + (mass2/(0.6*MSOLAR))**(1/3))
    I_orb = chirpMass**(5/3) / ((np.pi*freq0)**(4/3))
    #1PN freqD and freqDD
    fdot = fdot_pp * (1 + ((3*I_wd/I_orb)/(1 - (3*I_wd/I_orb))) )
    fddot = (11/3)*(fdot_pp**2/freq0 )* ((1 - (21/11)*(I_wd/I_orb)) / ((1 - (3*I_wd/I_orb))**3))
    fddot_0PN = (11/3)*(fdot ** 2) / freq0
    delta_fddot = fddot - fddot_0PN 
    fd_corr = fdot_pp * (3*I_wd/I_orb)/(1 - (3*I_wd/I_orb)) 

    #frequency bin
    df = (1/t_obs)
    dfd = (1/t_obs**2)
    dfdd = (1/ t_obs**3)

    """ print("TIDAL CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN*t_obs**3, "Hz")
    print("Tides Freq Derivative =", fdot, "Hz")
    print("Tides Freq Derivative Correction=", fd_corr * t_obs**2, "Hz")
    print("Tides Freq Second Derivative =", fddot, "Hz")
    print("The Tidal Fdd correction is:", delta_fddot * t_obs**3, "Hz")
    print("Tides Freq Derivative Dimensionless=", fdot * t_obs**2, "Hz")
    print("Tides Freq Second Derivative Dimensionless=", fddot * t_obs**3, "Hz")
    print("-----------------------------------------------")
    print("Frequency bin:", df)
    print("Change in Freq bin due to fdot:", dfd)
    print("Change in Freq bin due to fddot:", dfdd)
    print("-----------------------------------------------") """
    return fdot, fddot


def getFrequency_ChirpTotalMass(freq0, params, true_frequencies):
    chirpMass, totalMass = params
    # chirpMass = max(chirpMass, 1.e-10)
    # totalMass = max(totalMass, 1.e-10)

    #Scale values closer to unity - mHz/ ms
    # scale = 1.e3
    # chirpMass *= scale
    # totalMass *= scale
    # freq0 *= scale
    
    eta = ( chirpMass / totalMass )**(5/3)

    fdot_pp = (96/5)*(np.pi**(8/3))*(freq0**(11/3))*(chirpMass**(5/3))
    fdot = fdot_pp * (1 + ((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))

    fddot = fdot_pp * (fdot/freq0) * ((11/3) + (13/3)*
                                      ((743/1344)-(11*eta/16))*((8*np.pi*totalMass*freq0)**(2/3)))

    scale = 1.e0
    fdot *= scale
    fddot *= scale
    fdot_true = true_frequencies[0] * scale
    fddot_true = true_frequencies[1] * scale

    return [fdot - fdot_true, fddot - fddot_true]


def getRootFinder(freq0, fd, fdd, mass1, mass2):
    #root finding method for Mt and Mc from the 1PN GW equations

    mass1_exact = 0.7*MSOLAR
    mass2_exact = 0.6*MSOLAR

    chirpMass_guess = getChirpMass(mass1, mass2)
    totalMass_guess = getTotalMass(mass1, mass2)

    chirpMass_exact = getChirpMass(mass1_exact, mass2_exact)
    totalMass_exact = getTotalMass(mass1_exact, mass2_exact)

    params_guess = np.array([chirpMass_guess, totalMass_guess])
    exact = [chirpMass_exact, totalMass_exact]

    #do the iteration
    fx = lambda p : getFrequency_ChirpTotalMass(freq0, p, [fd, fdd])

    finalguess = optimize.fsolve(fx, params_guess)
    print(np.isclose([fx(finalguess)], [0.0, 0.0]))

    print("Final Guess:", finalguess)
    print("Real Values:", exact)
    return

def getNewtonRaphson(fx, fx_t, jacobian, guess, eps, maxiter):
    i = 0

    while True:
        i += 1
        jacobianX = jacobian(guess)
        hessianX = torch.autograd.functional.hessian(fx_t, (torch.tensor(guess[0]), torch.tensor(guess[1])), create_graph=True)

        delta = np.linalg.multi_dot([np.linalg.inv(hessianX), jacobianX])
        guess = guess - delta

        if np.linalg.norm(delta) < eps or i >= maxiter:
            if i >= maxiter:
                i = -1
            break

    return guess, i 