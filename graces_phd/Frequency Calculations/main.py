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

    params_true = np.array([getChirpMass(mass1, mass2), getTotalMass(mass1, mass2)])

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

    #Post Processing Function Calls
    #frequencyPostProcessing(freq0, t_obs, mass1, mass2)  
    #massFiltering("chirp mass (unstripped) 50000.txt")
    #makeHistogramPlots("total mass stripped 50000.txt", "Mt", 1, 2)
    #makeHistogramPlots("chirp mass stripped 50000.txt", "Mc", 0.55, 0.58)
    #makeCornerPlot("chirp mass stripped 50000.txt", 'total mass stripped 50000.txt', params_true)

    ########## TESTING CHECKS ###############
    chirpMass_test, totalMass_test = frequencyPostProcessing(freq0, t_obs, mass1, mass2)
    test_beta = []
    test_delta = []
    for i in range(len(chirpMass_test)):
        b, d, dd = Frequency_Tides(freq0, chirpMass_test[i]*MSOLAR, totalMass_test[i]*MSOLAR, t_obs)
        test_beta.append(b*(t_obs**2))
        test_delta.append(dd*(t_obs**3))
    #b, d, dd = Frequency_Tides(freq0, 0.5640277240412154*MSOLAR, 1.2943656723814734*MSOLAR, t_obs)
    print(test_beta)
    print(test_delta)

main(20.e-3, 0.7*MSOLAR, 0.6*MSOLAR, 7.6e-22*KPCSEC, 4.0*SECSYEAR)