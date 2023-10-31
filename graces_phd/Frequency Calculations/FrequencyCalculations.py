import math
from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg, optimize
import matplotlib.pyplot as plt
from HelperCalculations import *
import sympy as sp


def Frequency_1PN(freq0, mass1, mass2, t_obs):

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
    fddot_0PN = (11/3)*(fdot_pp**2) / freq0
    fddot_1PN = fddot_0PN*( (3/11)*(13/3)*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fdddot_1PN = (19/3) * ((fdot * fddot) / freq0) #* ( 1 + (2/19) * (fdot_pp / fdot) * (1 + ((13/3)*(fdot**2 / (freq0 * fddot)))) * (((743/1344)-(11*eta/16))*((8*np.pi*totalMass*freq0)**(2/3))) )

    delta_fddot_v1 = fddot - (11/3)* (fdot**2 / freq0)
    delta_fddot = fddot - fddot_0PN 

    #amp = np.pi**2/3 * chirpMass**(5/3) * freq0**(2/3) / dl

    #frequency bin
    df = (1/t_obs)
    dfd = (1/t_obs**2)
    dfdd = (1/ t_obs**3)

    # print("1PN CORRECTION TERMS")
    # print("-----------------------------------------------")
    # print("1PN Alpha =", freq0 * (t_obs))
    # print("1PN Beta =", fdot * (t_obs)**2)
    # print("1PN Gamma =", fddot * (t_obs)**3)
    # print("1PN delta:", delta_fddot_v1* (t_obs)**3)
    # print("1PN fdddot:", fdddot_1PN)
    # print("1PN fdddot Unitless:", fdddot_1PN * t_obs**4)
    return fdot, fddot, delta_fddot_v1 

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
    #Tides freqD and freqDD
    fdot = fdot_pp * (1 + ((3*I_wd/I_orb)/(1 - (3*I_wd/I_orb))) )
    fddot = (11/3)*(fdot_pp**2/freq0 )* (1 + (((26/11)*(3*I_wd/I_orb)) / (1 - (3*I_wd/I_orb))) + ( (19/11) * ((3*I_wd/I_orb) / (1 - (3*I_wd/I_orb)))**2 ))
    fdddot = (19/3) * ((fdot * fddot) / freq0) * (1 + (12/19)*((3*I_wd/I_orb)/(1 - (3*I_wd/I_orb)))*(1 - (7/9)*((fdot_pp**2)/fdot)*((1 - (3*I_wd/I_orb))**(-2)) ) )
    
    fddot_0PN = (11/3)*(fdot ** 2) / freq0
    fddot_0PN_v2 = (11/3)*(fdot_pp ** 2) / freq0
    delta_fddot = fddot - fddot_0PN 
    delta_fddot_v2 = fddot - fddot_0PN_v2
    fd_corr = fdot_pp * (3*I_wd/I_orb)/(1 - (3*I_wd/I_orb)) 

    #frequency bin
    df = (1/t_obs)
    dfd = (1/t_obs**2)
    dfdd = (1/ t_obs**3)

    #print("TIDAL CORRECTION TERMS")
    #print("-----------------------------------------------")
    # print("Tides Alpha =", freq0 * (t_obs))
    # print("Tides Beta =", fdot * (t_obs)**2)
    # print("Tides Gamma =", fddot * (t_obs)**3)
    # print("Tides Delta:", delta_fddot * (t_obs)**3)
    # print("moment of inertia:", I_wd)
    # print("Iorb:", I_orb)
    # print("Tides fdddot:", fdddot)
    # print("Tides Kappa:", fdddot * t_obs**4)

    return fdot, fddot, delta_fddot


def getFrequency_ChirpTotalMass_1PN(freq0, params, results_exact):
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

def getFrequency_ComponentMass(freq0, t_obs, params, results_exact):
    mass1, mass2 = params

    mass1 = max(mass1, 1.e-30)
    mass2 = max(mass2, 1.e-30)
    
    chirpMass = ((mass1*mass2)**(3/5)) / ((mass1 + mass2)**(1/5))

    fdot_pp = (96 * (np.pi**(8/3)) * (chirpMass**(5/3)) * (freq0**(11/3))) / 5
    
    I_wd = 8.51e-10 * ((mass1/(0.6*MSOLAR))**(1/3) + (mass2/(0.6*MSOLAR))**(1/3))
    I_orb = chirpMass**(5/3) / ((np.pi*freq0)**(4/3))


    fdot = fdot_pp * (1 + ((3*I_wd/I_orb)/(1 - (3*I_wd/I_orb))) )
    fddot = (11/3)*(fdot_pp**2/freq0 )* (1 + (((26/11)*(3*I_wd/I_orb)) / (1 - (3*I_wd/I_orb))) + ( (19/11) * ((3*I_wd/I_orb) / (1 - (3*I_wd/I_orb)))**2 ))

    #parameterize frequencies
    alpha = freq0*t_obs
    beta = fdot*(t_obs**2)
    gamma = fddot*(t_obs**3)
    delta = (gamma - (11/3)*((beta**2) / (alpha)))
    fdot_isNaN = np.isnan(fdot)
    fddot_isNaN = np.isnan(fddot)

    F=np.zeros(2)
    F[0] = beta - results_exact[0]
    F[1] = delta - results_exact[1]
    return F

def getFrequency_McMt_Tides(freq0, t_obs, params, results_exact):
    chirpMass, totalMass = params
    eta = ( (chirpMass/MSOLAR) / (totalMass/MSOLAR) )**(5/3)
    dm = (1-(4*eta))**(1/2)

    if math.isnan(dm):
        dm = 1.e-30
    
    mass1 = totalMass * (1 + dm) / 2
    mass2 = totalMass * (1 - dm) / 2

    fdot_pp = (96 * (np.pi**(8/3)) * (chirpMass**(5/3)) * (freq0**(11/3))) / 5
    
    I_wd = 8.51e-10 * ((mass1/(0.6*MSOLAR))**(1/3) + (mass2/(0.6*MSOLAR))**(1/3))
    I_orb = chirpMass**(5/3) / ((np.pi*freq0)**(4/3))


    fdot = fdot_pp * (1 + ((3*I_wd/I_orb)/(1 - (3*I_wd/I_orb))) )
    fddot = (11/3)*(fdot_pp**2/freq0 )* (1 + (((26/11)*(3*I_wd/I_orb)) / (1 - (3*I_wd/I_orb))) + ( (19/11) * ((3*I_wd/I_orb) / (1 - (3*I_wd/I_orb)))**2 ))

    #parameterize frequencies
    alpha = freq0*t_obs
    beta = fdot*(t_obs**2)
    gamma = fddot*(t_obs**3)
    delta = (gamma - (11/3)*((beta**2) / (alpha)))

    F=np.zeros(2)
    F[0] = (beta/ results_exact[0]) - 1
    F[1] = delta - results_exact[1]
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
                fx = lambda params: getFrequency_ChirpTotalMass_1PN(freq0, params, results_exact)[i]
                dF = lambda x, params: derivative(fx, params, step_size, x)
                jacobian[i, j] = dF(j, p)
        return jacobian

    #do the iteration
    fx = lambda p : getFrequency_ChirpTotalMass_1PN(freq0, p, results_exact)

    final_guess = optimize.newton(fx, params_guess, tol=1.e-15, maxiter=10000)

    # for i in range(len(final_guess)):
    #     final_guess[i] /= MSOLAR
    #     params_exact[i] /= MSOLAR

    print("Final Guess:", final_guess)
    print(fx(final_guess))
    print("Real Values:", params_exact )
    print(fx(params_exact))
    return

def getRootFinder_tides_componentMass(freq0, fdot, fddot, t_obs, mass1_exact, mass2_exact, mass1_guess, mass2_guess):
    #root finding method for Mt and Mc from the tide GW equations
    #beta = fdot*(t_obs**2)
    #delta = (fddot - (11/3)*(fdot**2/freq0))*t_obs**3

    params_guess = [mass1_guess, mass2_guess]

    params_exact = [mass1_exact, mass2_exact]

    results_exact = [fdot, fddot]

    #do the iteration
    fx = lambda p : getFrequency_ComponentMass(freq0, t_obs, p, results_exact)

    final_guess = optimize.newton(fx, params_guess, tol=1.e-15, maxiter=1000000)
    
    for i in range(len(final_guess)):
        final_guess[i] /= MSOLAR
        params_exact[i] /= MSOLAR

    final_guess_chirp = getChirpMass(final_guess[0], final_guess[1])
    final_guess_total = getTotalMass(final_guess[0], final_guess[1])

    # print("Final Guess (m1, m2):", final_guess)
    # print("Real Values (m2, m1):", params_exact )

    # print("Final Guess Chirp:", final_guess_chirp)
    # print("Real Chirp:", (getChirpMass(mass1_exact, mass2_exact))/MSOLAR)

    # print("Final Guess Total:", final_guess_total)
    # print("Real Total:", (getTotalMass(mass1_exact, mass2_exact))/MSOLAR)

    return final_guess_chirp, final_guess_total

def getRootFinder_tides_chirpTotalMass(freq0, beta, delta, t_obs, chirp_exact, total_exact, chirp_guess, total_guess):
    #root finding method for Mt and Mc from the tide GW equations
    params_guess = [chirp_guess, total_guess]

    params_exact = [chirp_exact, total_exact]

    results_exact = [beta, delta]

    #do the iteration
    fx = lambda p : getFrequency_McMt_Tides(freq0, t_obs, p, results_exact)

    final_guess = optimize.root(fx, params_guess, tol=1.e-10)
    
    for i in range(len(final_guess.x)):
        final_guess.x[i] /= MSOLAR
        params_exact[i] /= MSOLAR

    #print("Final Guess (Mc,Mt):", final_guess)
    #print("Real Values (Mc, Mt):", params_exact )

    return final_guess