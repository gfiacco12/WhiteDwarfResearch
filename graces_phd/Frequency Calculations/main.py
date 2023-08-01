import numpy as np
from Calculations import Frequency_1PN
from const import *

def main(freq0, chirpMass, totalMass, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    Frequency_1PN(freq0, chirpMass, totalMass, t_obs)

main(10.e-3, 0.522*MSOLAR, 1.4*MSOLAR, 4)