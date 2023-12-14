import numpy as np
from FrequencyCalculations import *
from HelperCalculations import calculateAmplitude, calculateAmplitude_phys, dataFFT
from const import *
from graphing import *
from postprocessing import *
from priors import resampling, drawSamples_M1M2, resampling_test

def main(freq0, mass1, mass2, dl, t_obs):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA
    #INPUTS: Initial frequency (Hz), chirp mass (s), total mass (s)

    #frequency calculations
    fD_1PN, fDD_1PN, delta= Frequency_1PN(freq0, mass1, mass2, t_obs)
    beta, delta = Frequency_Tides_Masses(freq0, [mass1, mass2], t_obs)

    # Some amplitude/SNR calculations
    amp = calculateAmplitude(1000, t_obs)
    amp_phys = calculateAmplitude_phys(dl, 0.564*MSOLAR, freq0)

    #makeDeltaPlot(freq0, mass1, mass2, t_obs, amp)

    #Post Processing Function Calls
    frequencyPostProcessing(freq0, "beta delta mass check 50000 v2.txt", t_obs, mass1, mass2, beta, delta)  

    #prior transfers
    #drawSamples_M1M2(10000, [0.3,0.3], [1.4,1.4], [mass1, mass2], freq0)
    #sigmas = [0.017, 0.0403]
    #resampling_test([beta, delta], freq0, t_obs, 500000, sigmas)
    #resampling([beta, delta], freq0, t_obs, [mass1, mass2])

main(20.e-3, 0.7*MSOLAR, 0.6*MSOLAR, 7.6e-22*KPCSEC, 4.0*SECSYEAR)