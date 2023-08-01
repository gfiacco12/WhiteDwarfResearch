from numpy import number
import numpy as np
from const import *

def Frequency_1PN(freq0, chirpMass, totalMass):
    eta = (chirpMass/totalMass)**(5/3)
    fdot_pp = 96/5*np.pi**(8/3)*freq0**(11/3)*chirpMass**(5/3)
    fdot = fdot_pp * (1 + ((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    fddot = fdot_pp*fdot * ((11/3) + (13/3)*((743/1344)-(11*eta/16))*(8*np.pi*totalMass*freq0)**(2/3))
    return fdot, fddot 