from numpy import number
import numpy as np
from const import *

def Frequency_1PN(freq0, chirpMass, totalMass, t_obs):

    eta = (chirpMass/totalMass)**(5/3)
    #0PN point particle freqD
    fdot_pp = 96/5*np.pi**(8/3)*freq0**(11/3)*chirpMass**(5/3)
    
    #1PN freqD and freqDD
    fdot = fdot_pp * (1 + ((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fddot = fdot_pp * (fdot/freq0) * ((11/3) + (13/3)*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fddot_0PN = (fdot ** 2) / freq0
    delta_fddot = fddot - fddot_0PN 

    #frequency bin
    df = (1/t_obs)
    df_1 = 3.1e-8 * (1/t_obs)

    print("0PN point particle fdot:", fdot_pp, "Hz")
    print("0PN freqDD:", fddot_0PN, "Hz")
    print("1PN Freq Derivative =", fdot, "Hz")
    print("1PN Freq Second Derivative =", fddot, "Hz")
    print("The 1PN Fdd correction is:", delta_fddot, "Hz")
    print("Fractional 1PN Fdd:", 1 - delta_fddot/fddot, "Hz")
    print("-----------------------------------------------")
    print("Frequency bin:", df_1, "Hz")
    return fdot, fddot, delta_fddot 