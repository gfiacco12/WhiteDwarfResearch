#post-processing code for samples from mcmc
from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg, optimize
import matplotlib.pyplot as plt
from HelperCalculations import *
import corner
from FrequencyCalculations import *
import math


def frequencyPostProcessing(freq0, t_obs, mass1, mass2):
#start by reading in txt file
    beta = []
    delta = []

    f = open("mcmc samples (alpha beta delta) 50000 steps.txt", 'r')
    
    for row in f:
        row = row.split()
        for i in range(len(row)): 
            if i == 2:
                beta.append(float(row[i]))
            if i == 3:
                delta.append(float(row[i]))

    #now convert to masses using root finder
    chirpMass = []
    totalMass = []
    for i in range(len(beta)):
        chirpMass_guess, totalMass_guess = getRootFinder_tides(freq0, beta[i], delta[i], t_obs, mass1, mass2, 0.6*MSOLAR, 0.6*MSOLAR)
        chirpMass.append(chirpMass_guess)
        totalMass.append(totalMass_guess)

    #remove nans and negative values
    chirpMass_new = []
    totalMass_new = []

    for mass in range(len(chirpMass)):
        if not math.isnan(chirpMass[mass]) and totalMass[mass] > 0:
            chirpMass_new.append(chirpMass[mass])
            totalMass_new.append(totalMass[mass])          
    
    print("Final Guess Chirp:", np.mean(chirpMass_new))
    realChirp = (getChirpMass(mass1, mass2))/MSOLAR
    print("Real Chirp:", realChirp)

    print("Final Guess Total:", np.mean(totalMass_new))
    realTotal = (getTotalMass(mass1, mass2))/MSOLAR
    print("Real Total:", realTotal)
    
    