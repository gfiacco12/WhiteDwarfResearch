from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg, optimize
import matplotlib.pyplot as plt
from HelperCalculations import *
import sympy as sp


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
    fddot = fdot_pp * (fdot/freq0) * ((11/3) + (13/3)*((743/1344)-(11*eta/16))*((8*np.pi*totalMass*freq0)**(2/3)))
    fddot_0PN = (11/3)*(fdot * fdot_pp) / freq0
    fddot_1PN = fddot_0PN*( (3/11)*(13/3)*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    
    
    delta_fddot = fddot - fddot_0PN 

    amp = np.pi**2/3 * chirpMass**(5/3) * freq0**(2/3) / dl

    #frequency bin
    df = (1/t_obs)
    dfd = (1/t_obs**2)
    dfdd = (1/ t_obs**3)

    # print("1PN CORRECTION TERMS")
    # print("-----------------------------------------------")
    # print("eta:", eta)
    # print("0PN point particle fdot:", fdot_pp * t_obs**2, "Hz")
    # print("1PN Freq Derivative =", fdot, "Hz")
    # print("1PN Freq Second Derivative =", fddot, "Hz")
    # print("Test 1PN fdd corr:", fddot_1PN)
    # print("-----------------------------------------------")
    # print("Frequency bin:", df)
    # print("Change in Freq bin due to fdot:", dfd)
    # print("Change in Freq bin due to fddot:", dfdd)
    # print("-----------------------------------------------")
    return fdot, fddot

def Frequency_Tides(freq0, chirpMass, totalMass, t_obs):
    eta = (chirpMass/totalMass)**(5/3)
    deltaM = np.sqrt(1 - 4*eta)
    mass1 = totalMass * (1 + deltaM) / 2
    mass2 = totalMass * (1 - deltaM) / 2
    return Frequency_Tides_Internal(freq0, chirpMass, mass1, mass2, t_obs)

def Frequency_Tides_Masses(freq0, mass1, mass2, t_obs):
    chirpMass = getChirpMass(mass1, mass2)
    return Frequency_Tides_Internal(freq0, chirpMass, mass1, mass2, t_obs)
    
def Frequency_Tides_Internal(freq0,chirpMass, mass1, mass2, t_obs):
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

    #print("TIDAL CORRECTION TERMS")
    #print("-----------------------------------------------")
    #print("eta:", eta)
    #print("0PN point particle fdot:", fdot_pp, "Hz")
    #print("0PN freqDD:", fddot_0PN*t_obs**3, "Hz")
    #print("Tides Freq Derivative =", fdot, "Hz")
    #print("Tides Freq Derivative Correction=", fd_corr * t_obs**2, "Hz")
    #print("Tides Freq Second Derivative =", fddot, "Hz")
    #print("The Tidal Fdd correction is:", delta_fddot * t_obs**3, "Hz")
    #print("Tides Freq Derivative Dimensionless=", fdot * t_obs**2, "Hz")
    #print("Tides Freq Second Derivative Dimensionless=", fddot * t_obs**3, "Hz")
    #print("-----------------------------------------------")
    #print("Frequency bin:", df)
    #print("Change in Freq bin due to fdot:", dfd)
    #print("Change in Freq bin due to fddot:", dfdd)
    #print("-----------------------------------------------") 
    return fdot, fddot


def getFrequency_ChirpTotalMass(freq0, params, results_exact):
    chirpMass, totalMass = params
    
    eta = ( chirpMass / totalMass )**(5/3)

    fdot_pp = (96 * (np.pi**(8/3)) * (chirpMass**(5/3)) * (freq0**(11/3))) / 5

    correction_1PN = ((743/1344) - ((11 * eta) / 16)) * ((8*np.pi*totalMass*freq0)**(2/3))

    fdot = fdot_pp * (1 + correction_1PN)
    fddot = ((fdot_pp * fdot) / freq0) * ((11/3) + ((13/3) * correction_1PN))

    F=np.zeros(2)
    F[0] = fdot - results_exact[0]
    F[1] = fddot - results_exact[1]
    return F

def getFrequency_ComponentMass(freq0, params, results_exact):
    mass1, mass2 = params

    mass1 = max(mass1, 1.e-30)
    mass2 = max(mass2, 1.e-30)
    
    chirpMass = ((mass1*mass2)**(3/5)) / ((mass1 + mass2)**(1/5))

    fdot_pp = (96 * (np.pi**(8/3)) * (chirpMass**(5/3)) * (freq0**(11/3))) / 5
    
    I_wd = 8.51e-10 * ((mass1/(0.6*MSOLAR))**(1/3) + (mass2/(0.6*MSOLAR))**(1/3))
    I_orb = chirpMass**(5/3) / ((np.pi*freq0)**(4/3))


    fdot = fdot_pp * (1 + ((3*I_wd/I_orb)/(1 - (3*I_wd/I_orb))) )
    fddot = (11/3)*(fdot_pp**2/freq0 )* ((1 - (21/11)*(I_wd/I_orb)) / ((1 - (3*I_wd/I_orb))**3))

    fdot_isNaN = np.isnan(fdot)
    fddot_isNaN = np.isnan(fddot)

    F=np.zeros(2)
    F[0] = fdot - results_exact[0]
    F[1] = fddot - results_exact[1]
    return F

def getRootFinder_1PN(freq0, fdot, fddot, mass1_exact, mass2_exact, mass1_guess, mass2_guess):
    #root finding method for Mt and Mc from the 1PN GW equations

    chirpMass_guess = getChirpMass(mass1_guess, mass2_guess)
    totalMass_guess = getTotalMass(mass1_guess, mass2_guess)
    params_guess = [chirpMass_guess, totalMass_guess]

    chirpMass_exact = getChirpMass(mass1_exact, mass2_exact)
    totalMass_exact = getTotalMass(mass1_exact, mass2_exact)
    params_exact = [chirpMass_exact, totalMass_exact]

    results_exact = [fdot, fddot]

    step_size = [1.e-8, 1.e-8]
    def get_Jacobian(p):
        jacobian = np.zeros((np.size(p), np.size(p)))
        for i in range(len(p)):
            for j in range(len(p)):
                fx = lambda params: getFrequency_ChirpTotalMass(freq0, params, results_exact)[i]
                dF = lambda x, params: derivative(fx, params, step_size, x)
                jacobian[i, j] = dF(j, p)
        return jacobian

    #do the iteration
    fx = lambda p : getFrequency_ChirpTotalMass(freq0, p, results_exact)

    final_guess = optimize.newton(fx, params_guess, tol=1.e-15, maxiter=10000)

    # for i in range(len(final_guess)):
    #     final_guess[i] /= MSOLAR
    #     params_exact[i] /= MSOLAR

    print("Final Guess:", final_guess)
    print(fx(final_guess))
    print("Real Values:", params_exact )
    print(fx(params_exact))
    return

def getRootFinder_tides(freq0, fdot, fddot, mass1_exact, mass2_exact, mass1_guess, mass2_guess):
    #root finding method for Mt and Mc from the 1PN GW equations

    params_guess = [mass1_guess, mass2_guess]

    params_exact = [mass1_exact, mass2_exact]

    results_exact = [fdot, fddot]

    #do the iteration
    fx = lambda p : getFrequency_ComponentMass(freq0, p, results_exact)

    final_guess = optimize.newton(fx, params_guess, tol=1.e-15, maxiter=1000000)

    for i in range(len(final_guess)):
        final_guess[i] /= MSOLAR
        params_exact[i] /= MSOLAR

    print("Final Guess:", final_guess)
    # print(fx(final_guess))
    print("Real Values:", params_exact )
    # print(fx(params_exact))
    return