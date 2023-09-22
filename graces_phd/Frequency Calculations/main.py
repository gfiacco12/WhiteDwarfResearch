import numpy as np
from FrequencyCalculations import Frequency_1PN, Frequency_Tides, Frequency_Tides_Masses, getRootFinder_1PN, getRootFinder_tides
from FisherCalculations import *
from HelperCalculations import calculateAmplitude, calculateAmplitude_phys, getChirpMass, getTotalMass
from const import *
from graphing import makeDeltaPlot

def main(freq0, mass1, mass2, dl, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    #frequency calculations
    fD_1PN, fDD_1PN, deltafDD_1PN = Frequency_1PN(freq0, mass1, mass2, t_obs)
    fD_tide, fDD_tide, deltafDD_tide = Frequency_Tides_Masses(freq0, mass1, mass2, t_obs)
    
    #getRootFinder_1PN(freq0, fD_1PN, fDD_1PN, mass1, mass2, 0.6*MSOLAR, 0.6*MSOLAR)
    #getRootFinder_tides(freq0, fD_tide, fDD_tide, mass1, mass2, 0.6*MSOLAR, 0.6*MSOLAR)

    # Some amplitude/SNR calculations
    amp = calculateAmplitude(1000, t_obs)
    amp_phys = calculateAmplitude_phys(dl, 0.564*MSOLAR, freq0)
    snr = getSNR(t_obs, amp_phys, 1.5, freq0, fD_1PN, fDD_1PN)

    #fisher_frequencies(freq0, fD_1PN, fDD_1PN, t_obs, amp)
    #fisher_frequencies(freq0, fD_tide, fDD_tide, t_obs, amp)
    #fisher_chirpTotal(freq0, mass1, mass2, t_obs, amp)
    #fisher_masses(freq0, mass1, mass2, t_obs, amp)

    makeDeltaPlot(freq0, mass1, mass2, t_obs, amp)
main(20.e-3, 0.7*MSOLAR, 0.6*MSOLAR, 7.6e-22*KPCSEC, 4*SECSYEAR)