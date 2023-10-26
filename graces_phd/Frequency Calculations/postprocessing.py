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
    params_true = np.array([getChirpMass(mass1, mass2)/MSOLAR, getTotalMass(mass1, mass2)/MSOLAR])
    for i in range(len(beta)):
        chirpMass_guess, totalMass_guess = getRootFinder_tides_chirpTotalMass(freq0, beta[i], delta[i], t_obs, params_true[0], params_true[1], 0.5*MSOLAR, 1.2*MSOLAR)
        chirpMass.append(chirpMass_guess)
        totalMass.append(totalMass_guess)
    np.savetxt('chirp & total masses not stripped.txt', np.array([chirpMass, totalMass]))
    #remove nans and negative values
    chirpMass_new = []
    totalMass_new = []

    # for mass in range(len(chirpMass)):
    #     if not math.isnan(chirpMass[mass]) and totalMass[mass] > 0:
    #         chirpMass_new.append(chirpMass[mass])
    #         totalMass_new.append(totalMass[mass])          
    
    print("Final Guess Chirp:", np.mean(chirpMass_new))
    realChirp = (getChirpMass(mass1, mass2))/MSOLAR
    print("Real Chirp:", realChirp)

    print("Final Guess Total:", np.mean(totalMass_new))
    realTotal = (getTotalMass(mass1, mass2))/MSOLAR
    print("Real Total:", realTotal)
    
def massFiltering(data_file):
    
    masses = []
    mass_new = []

    f = open(data_file, 'r')

    for row in f:
        elements = row.split(' ')
        elements = list(map(lambda e : float(e), elements))
        masses += elements

    for mass in range(len(masses)):
        if masses[mass] > 0:
            mass_new.append(masses[mass])
    print(len(masses))
    print(len(mass_new))
    np.savetxt("total mass filtered (only neg).txt", mass_new)