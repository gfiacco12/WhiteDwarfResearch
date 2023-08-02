import numpy as np
from Calculations import Frequency_1PN, Frequency_Tides
from const import *

def main(freq0, mass1, mass2, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    Frequency_1PN(freq0, mass1, mass2, t_obs)
    Frequency_Tides(freq0, mass1, mass2, t_obs)

main(10.e-3, 0.6*MSOLAR, 1.0*MSOLAR, 4)