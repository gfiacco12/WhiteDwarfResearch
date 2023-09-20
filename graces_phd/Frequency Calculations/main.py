import numpy as np
from FrequencyCalculations import Frequency_1PN, Frequency_Tides, Frequency_Tides_Masses, getRootFinder_1PN, getRootFinder_tides
from FisherCalculations import GWsignal_chirpTotal, GWsignal_masses, getFisherMatrix, getSNR, GWsignal
from HelperCalculations import calculateAmplitude, calculateAmplitude_phys, getChirpMass, getTotalMass
from const import *

def main(freq0, mass1, mass2, dl, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    #frequency calculations
    fD_1PN, fDD_1PN, deltafDD_1PN = Frequency_1PN(freq0, mass1, mass2, dl, t_obs)
    fD_tide, fDD_tide = Frequency_Tides_Masses(freq0, mass1, mass2, t_obs)
    
    #getRootFinder_1PN(freq0, fD_1PN, fDD_1PN, mass1, mass2, 0.6*MSOLAR, 0.6*MSOLAR)
    #getRootFinder_tides(freq0, fD_tide, fDD_tide, mass1, mass2, 0.6*MSOLAR, 0.6*MSOLAR)

    # Some amplitude/SNR calculations
    amp = calculateAmplitude(100, t_obs)
    amp_phys = calculateAmplitude_phys(dl, 0.564*MSOLAR, freq0)
    snr = getSNR(t_obs, amp_phys, 1.5, freq0, fD_1PN, fDD_1PN)

    def fisher_frequencies(fd, fdd):
        fx = lambda t, p : GWsignal(t, t_obs, p)
        alpha = freq0 * t_obs
        beta = fd * (t_obs**2)
        gamma = fdd * (t_obs**3)
        
        getFisherMatrix(
            t_obs, 
            fx, 
            np.array([amp, 1.5, alpha, beta, gamma]),
            [r'$A$', r'$\phi_{0}$',  r'$\alpha}$', r'$\beta$', r'$\gamma$']
        )
        delta = (gamma - (11 / 3)*(beta**2 / alpha))
        print("delta:", delta)
    #fisher_frequencies(fD_1PN, fDD_1PN)
    #fisher_frequencies(fD_tide, fDD_tide)

    def fisher_chirpTotal():
        fx = lambda t, p : GWsignal_chirpTotal(t, t_obs, p)
        chirpMass = getChirpMass(mass1, mass2)
        totalMass = getTotalMass(mass1, mass2)
        getFisherMatrix(
            t_obs, 
            fx, 
            np.array([amp, 1.5, freq0, chirpMass, totalMass]),
            [r'$A$', r'$\phi_{0}$',  r'$f0$', r'$MC$', r'$MT$']
        )
    #fisher_chirpTotal()

    def fisher_masses():
        fx = lambda t, p : GWsignal_masses(t, t_obs, p)
        getFisherMatrix(
            t_obs, 
            fx, 
            np.array([amp, 1.5, freq0, mass1, mass2]),
            [r'$A$', r'$\phi_{0}$',  r'$f0$', r'$M1$', r'$M2$']
        )
    #fisher_masses()

main(20.e-3, 0.7*MSOLAR, 0.6*MSOLAR, 9e-20*KPCSEC, 4*SECSYEAR)