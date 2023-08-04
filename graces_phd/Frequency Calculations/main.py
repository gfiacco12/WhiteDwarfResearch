import numpy as np
from Calculations import Frequency_1PN, Frequency_Tides, getFisherMatrix, calculateAmplitude
from const import *

def main(freq0, mass1, mass2, dl, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    fD_1PN, fDD_1PN = Frequency_1PN(freq0, mass1, mass2, dl, t_obs)
    fD_tide, fDD_tide = Frequency_Tides(freq0, mass1, mass2, dl, t_obs)

    fish = getFisherMatrix(t_obs, calculateAmplitude(10, t_obs), 1.5, freq0, fD_1PN, fDD_1PN)
    print(fish)
main(10.e-3, 0.6*MSOLAR, 1.0*MSOLAR, 10*KPCSEC, 4*SECSYEAR)