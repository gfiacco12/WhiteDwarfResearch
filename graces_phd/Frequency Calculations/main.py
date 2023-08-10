import numpy as np
from FrequencyCalculations import Frequency_1PN, Frequency_Tides, getFxAndJacobian
from FisherCalculations import getFisherMatrix, fisher, getSNR
from HelperCalculations import calculateAmplitude, calculateAmplitude_phys
from const import *

def main(freq0, mass1, mass2, dl, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    #frequency calculations
    fD_1PN, fDD_1PN = Frequency_1PN(freq0, mass1, mass2, dl, t_obs)
    fD_tide, fDD_tide = Frequency_Tides(freq0, mass1, mass2, dl, t_obs)
    
    getFxAndJacobian(freq0, 0.6*MSOLAR, 0.6*MSOLAR)

    # Some amplitude/SNR calculations
    amp = calculateAmplitude(100, t_obs)
    amp_phys = calculateAmplitude_phys(dl, 0.522*MSOLAR, freq0)
    snr = getSNR(t_obs, amp_phys, 1.5, freq0, fD_1PN, fDD_1PN)

    getFisherMatrix(t_obs, amp, 1.5, freq0, fD_1PN, fDD_1PN)
    fisher(t_obs, amp, 1.5, freq0, fD_1PN, fDD_1PN, 1)

main(10.e-3, 0.6*MSOLAR, 0.7*MSOLAR, 5.6e-20*KPCSEC, 4*SECSYEAR)