from numpy import number
import numpy as np
from const import *

def Frequency_1PN(freq0, mass1, mass2, t_obs):

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

    #frequency bin
    df = (1/t_obs)
    df_1 = 3.1e-8 * (1/t_obs)

    print("1PN CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("1PN Freq Derivative =", fdot, "Hz")
    print("1PN Freq Second Derivative =", fddot, "Hz")
    print("1PN fdot corr:", fd_corr)
    print("The 1PN Fdd correction is:", delta_fddot, "Hz")
    print("Frequency bin:", df_1, "Hz")
    print("-----------------------------------------------")
    return fdot, fddot, delta_fddot 


def Frequency_Tides(freq0, mass1, mass2, t_obs):

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

    #frequency bin
    df = (1/t_obs)
    df_1 = 3.1e-8 * (1/t_obs)

    print("TIDAL CORRECTION TERMS")
    print("-----------------------------------------------")
    print("eta:", eta)
    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("Tides Freq Derivative =", fdot, "Hz")
    print("Tides Freq Second Derivative =", fddot, "Hz")
    print("Tidal fdot corr:", fd_corr)
    print("The Tidal Fdd correction is:", delta_fddot, "Hz")
    print("Frequency bin:", df_1, "Hz")
    print("-----------------------------------------------")
    return fdot, fddot, delta_fddot 


def fisherMatrix(t_obs):
    #define f(t), phi(t)
    #declare variables for fisher
    t = np.linspace[0, t_obs, 1000]

    #create a function to calculate h(t). Use LISA code to grab values for the params. Copy fisher matrix
    #code from my github to do the actual derivative and looping through the values and elements of matrix
    # maybe look into scipy function to take derivatives instead of coding it myself

    phi = lambda phi0, f0, fD, fDD : phi0 + 2*np.pi*f0*t + np.pi*fD*t**2 + (np.pi/3)*fDD*t**3
    h = lambda A, phi0, f0, fD, fDD : A * np.cos(2 * np.pi * phi(phi0, f0, fD, fDD))

    #take derivatives for elements of fisher matrix
    return

def fisherMatrixElement(noise = 1):
    const = 2 / noise
    
    return

def h(t):
    return

def add(x, y):
    return x, y

add = lambda x, y : x + y