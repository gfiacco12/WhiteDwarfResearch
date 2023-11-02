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
from graphing import *
import math


def frequencyPostProcessing(freq0, file_name, t_obs, mass1, mass2):
#start by reading in txt file
    beta = []
    delta = []

    f = open(file_name, 'r')
    
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
        #starting guess for chirp mass and total mass
        fdot = beta[i]/(t_obs**2)
        chirp_guess = (fdot * 5 / (96 * (np.pi**(8/3)) * (freq0**(11/3))))**(3./5)
        total_guess = chirp_guess / (0.24)**(3./5)
        final_guess = getRootFinder_tides_chirpTotalMass(freq0, beta[i], delta[i], t_obs, params_true[0], params_true[1], chirp_guess, total_guess)
        #filter the masses
        chirp = final_guess.x[0]
        total = final_guess.x[1]
        if final_guess.success == True:
            chirpMass.append(chirp)
            totalMass.append(total)
    ###########################################################
    print(len(chirpMass))
    print("Final Guess Chirp:", np.mean(chirpMass))
    realChirp = params_true[0]
    print("Real Chirp:", realChirp)

    print(len(totalMass))
    print("Final Guess Total:", np.mean(totalMass))
    realTotal = params_true[1]
    print("Real Total:", realTotal)
    
    #save new masses to file for plotting
    np.savetxt("chirp mass stripped 150000.txt", chirpMass)
    np.savetxt("total mass stripped 150000.txt", totalMass) 

    makeHistogramPlots(totalMass, "Mt (150000 steps)", 1, 2)
    makeHistogramPlots(chirpMass, "Mc (150000 steps)", 0.55, 0.58)
    makeCornerPlot(chirpMass, totalMass, params_true)
    return chirpMass, totalMass

def massFiltering(data_file, data):
    
    masses = []
    mass_new = []
    bad_mass = []

    f = open(data_file, 'r')

    for row in f:
        elements = row.split(' ')
        elements = list(map(lambda e : float(e), elements))
        masses += elements

    for mass in range(len(masses)):
        if data.success == True:
            mass_new.append(masses[mass])
        else:
            bad_mass.append(masses[mass])
    print(bad_mass)
    return mass_new