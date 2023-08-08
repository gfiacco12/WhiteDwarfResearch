from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg
import matplotlib.pyplot as plt
from HelperCalculations import *
from FisherCalculations import getParamsWithStep

def Frequency_1PN(freq0, mass1, mass2, dl, t_obs):

    chirpMass = getChirpMass(mass1, mass2)
    totalMass = getTotalMass(mass1, mass2)
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
    dfd = (1/t_obs**2)
    dfdd = (1/ t_obs**3)


    print("1PN CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("1PN Freq Derivative =", fdot, "Hz")
    print("1PN Freq Second Derivative =", fddot, "Hz")
    print("1PN fdot corr:", fd_corr)
    print("The 1PN Fdd correction is:", delta_fddot, "Hz")
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
    fddot = (11/3)*(fdot_pp * fdot/freq0 )* ((1 - (21/11)*(I_wd/I_orb)) / ((1 - (3*I_wd/I_orb))**2))
    fddot_0PN = (11/3)*(fdot ** 2) / freq0
    delta_fddot = fddot - fddot_0PN 
    fd_corr = (3*I_wd/I_orb)/(1 - (3*I_wd/I_orb)) 

    #frequency bin
    df = (1/t_obs)
    dfd = (1/t_obs**2)
    dfdd = (1/ t_obs**3)

    print("TIDAL CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("Tides Freq Derivative =", fdot, "Hz")
    print("Tides Freq Second Derivative =", fddot, "Hz")
    print("Tidal fdot corr:", fd_corr)
    print("The Tidal Fdd correction is:", delta_fddot, "Hz")
    print("-----------------------------------------------")
    print("Frequency bin:", df)
    print("Change in Freq bin due to fdot:", dfd)
    print("Change in Freq bin due to fddot:", dfdd)
    print("-----------------------------------------------")
    return fdot, fddot


def Newton_Raphson(freq0, mass1_theory, mass2_theory, mass1, mass2, dl, t_obs, eps):
    #root finding method for Mt and Mc from the 1PN GW equations
    #need initial freq, fD, fDD data, and your symmetric mass ratio guess (for now assume equal mass: eta = 1/4)
    #Need to define initial guesses first

    F = np.array([Frequency_1PN(freq0, mass1_theory, mass2_theory, dl, t_obs)[0], 
                  Frequency_1PN(freq0, mass1_theory, mass2_theory, dl, t_obs)[1]])
    
    #Jacobian - issue, freq functions take m1 m2 as inputs, I need a list of params as inputs to take derivatives
    #similar to fisher matrix calculation

    #Real param values
    chirpMass = getChirpMass(mass1, mass2)
    totalMass = getTotalMass(mass1, mass2)

    params = [chirpMass, totalMass]
    
    #starting guess param values
    eta = (mass1_theory * mass2_theory) / (mass1_theory + mass2_theory)**2
    chirpMass_guess = ((5 / (96 * np.pi**(8/3))) * (Frequency_1PN(freq0, mass1_theory, mass2_theory, dl, t_obs)[0] / freq0**(11/3)))**(3/5)
    totalMass_guess = chirpMass_guess / eta**(3/5)

    F_norm = np.linalg.norm(F, ord=2)


    return