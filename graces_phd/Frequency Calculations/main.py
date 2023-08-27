import numpy as np
from FrequencyCalculations import Frequency_1PN, Frequency_Tides, getRootFinder
from FisherCalculations import getFisherMatrix, getSNR, GWsignal
from HelperCalculations import calculateAmplitude, calculateAmplitude_phys, getChirpMass
from const import *

def main(freq0, mass1, mass2, dl, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    #frequency calculations
    fD_1PN, fDD_1PN = Frequency_1PN(freq0, mass1, mass2, dl, t_obs)
    fD_tide, fDD_tide = Frequency_Tides(freq0, mass1, mass2, dl, t_obs)
    
    getRootFinder(freq0, fD_1PN, fDD_1PN, mass1, mass2, 0.6*MSOLAR, 0.6*MSOLAR)

    # Some amplitude/SNR calculations
    amp = calculateAmplitude(100, t_obs)
    amp_phys = calculateAmplitude_phys(dl, 0.564*MSOLAR, freq0)
    snr = getSNR(t_obs, amp_phys, 1.5, freq0, fD_1PN, fDD_1PN)
    
    """ print("Amplitude", amp)
    print("Amp phys:", amp_phys)
    print("SNR:", snr) 

    params = np.array([amp, 1.5, freq0, fD_1PN, fDD_1PN])
    t = t_obs/2
    h= GWsignal(t, t_obs, params) """
    
    #print("h(t)=", h)

    #getFisherMatrix(t_obs, amp, 1.5, freq0, fD_tide, fDD_tide)
    # print("delta alpha = ", (3**(1/2))/(8**(1/2)*np.pi*100))
    # print("delta beta = ", (5**(1/2))/(np.pi*100))
    # print("delta gamma = ", (3*(7**(1/2)))/(np.pi*100))

main(10.e-3, 0.7*MSOLAR, 0.6*MSOLAR, 9e-20*KPCSEC, 4*SECSYEAR)