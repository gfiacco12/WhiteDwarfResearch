import numpy as np
from FrequencyCalculations import *
from HelperCalculations import calculateAmplitude, calculateAmplitude_phys, getChirpMass, getTotalMass
from const import *
from graphing import *
from postprocessing import *

def main(freq0, mass1, mass2, dl, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    #frequency calculations
    fD_1PN, fDD_1PN, deltafDD_1PN = Frequency_1PN(freq0, mass1, mass2, t_obs)
    fD_tide, fDD_tide, deltafDD_tide = Frequency_Tides_Masses(freq0, mass1, mass2, t_obs)
    
    beta = fD_tide * (t_obs**2)
    delta = (fDD_tide - (11/3)*(fD_tide**2 / freq0))*(t_obs**3)

    params_true = np.array([getChirpMass(mass1, mass2), getTotalMass(mass1, mass2)])
    print(params_true/MSOLAR)
    #getRootFinder_tides_chirpTotalMass(freq0, beta, delta, t_obs, params_true[0], params_true[1], 0.3*MSOLAR, 1.4*MSOLAR)
    
    # Some amplitude/SNR calculations
    amp = calculateAmplitude(1000, t_obs)
    amp_phys = calculateAmplitude_phys(dl, 0.564*MSOLAR, freq0)
    snr = getSNR(t_obs, amp_phys, 1.5, freq0, fD_1PN, fDD_1PN)
    
    #fisher_frequencies(freq0, fD_1PN, fDD_1PN, t_obs, amp)
    #fisher_frequencies(freq0, fD_tide, fDD_tide, t_obs, amp)
    #fisher_chirpTotal(freq0, mass1, mass2, t_obs, amp)
    #fisher_masses(freq0, mass1, mass2, t_obs, amp)

    #makeDeltaPlot(freq0, mass1, mass2, t_obs, amp)

    frequencyPostProcessing(freq0, t_obs, mass1, mass2)
    # makeCornerPlots("chirp mass.txt", "Mc")
    # makeCornerPlots("chirp & total masses.txt", "Mt")
    
    #massFiltering("total mass not stripped.txt")
    # makeHistogramPlots("total masses.txt", "Mt")
    # makeCornerPlot("chirp mass.txt", 'total masses.txt', params_true)
main(20.e-3, 0.7*MSOLAR, 0.6*MSOLAR, 7.6e-22*KPCSEC, 4.0*SECSYEAR)