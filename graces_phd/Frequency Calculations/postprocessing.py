#post-processing code for samples from mcmc
import random
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
    ##### TESTING #####
    test_beta = []
    test_delta = []
    i = 0
    while i < 10:
        element = random.randrange(0, len(beta))
        test_beta.append(beta[element])
        test_delta.append(delta[element])
        i += 1
    print(test_beta)
    print(test_delta)
    #now convert to masses using root finder
    chirpMass = []
    totalMass = []
    params_true = np.array([getChirpMass(mass1, mass2)/MSOLAR, getTotalMass(mass1, mass2)/MSOLAR])
    for i in range(len(test_beta)):
        chirpMass_guess, totalMass_guess = getRootFinder_tides_chirpTotalMass(freq0, test_beta[i], test_delta[i], t_obs, params_true[0], params_true[1], 0.5*MSOLAR, 1.2*MSOLAR)
        chirpMass.append(chirpMass_guess)
        totalMass.append(totalMass_guess)
    print(chirpMass)
    print(totalMass)
    #np.savetxt('chirp mass (unstripped) 50000.txt', chirpMass)      
    #np.savetxt("total mass (unstripped) 50000.txt", totalMass)

    '''
    #remove nans and negatives
    newChirpMass = massFiltering("chirp mass (unstripped) 50000.txt")
    newTotalMass = massFiltering("total mass (unstripped) 50000.txt")

    print("Final Guess Chirp:", np.mean(newChirpMass))
    realChirp = params_true[0]
    print("Real Chirp:", realChirp)

    print("Final Guess Total:", np.mean(newTotalMass))
    realTotal = params_true[1]
    print("Real Total:", realTotal)
    
    #save new masses to file for plotting
    np.savetxt("chirp mass stripped 50000.txt", newChirpMass)
    np.savetxt("total mass stripped 50000.txt", newTotalMass) '''
    return chirpMass, totalMass

def massFiltering(data_file):
    
    masses = []
    mass_new = []

    f = open(data_file, 'r')

    for row in f:
        elements = row.split(' ')
        elements = list(map(lambda e : float(e), elements))
        masses += elements

    for mass in range(len(masses)):
        if masses[mass] > 0 and not math.isnan(masses[mass]):
            mass_new.append(masses[mass])
    return mass_new