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
    print("Final Guess Chirp:", np.mean(chirpMass))
    realChirp = params_true[0]
    print("Real Chirp:", realChirp)

    print("Final Guess Total:", np.mean(totalMass))
    realTotal = params_true[1]
    print("Real Total:", realTotal)

    makeCornerPlot(chirpMass, totalMass, params_true)

    postProcessMcMt_to_Components(chirpMass, totalMass, mass1, mass2)

    return chirpMass, totalMass


def postProcessMcMt_to_Components(chirpmass, totalmass, mass1, mass2):
    convertedM1 = []
    convertedM2 = []
    for i in range(len(chirpmass)):
        m1, m2 = get_comp_mass_Mc_Mt(chirpmass[i], totalmass[i])
        if m1 < 1.4 and m2 < 1.4:
            if np.isnan(m1) == False and np.isnan(m2) == False:
                convertedM1.append(m1)
                convertedM2.append(m2)

    print("Final M1:", np.average(convertedM1))
    print("Final M2:", np.average(convertedM2))

    bloop = []
    for i in range(len(convertedM1)):
        bloop.append(np.average([convertedM1[i], convertedM2[i]]))

    makeHistogramPlots([x - 0.7 for x in convertedM1], "M1 (150000 steps)", -1, 1)
    makeHistogramPlots([x - 0.6 for x in convertedM2], "M2 (150000 steps)", -1, 1)
    makeHistogramPlots(bloop, "bloopity bloop bloop (150000 steps)", 0.6, 1)
    makeCornerPlot(convertedM1, convertedM2, np.array([mass1/MSOLAR, mass2/MSOLAR]))
    return
    

def betaDeltaM1M2Converter(beta, delta, freq0, t_obs, m1, m2):
    #now convert to masses using root finder
    chirpMass = []
    totalMass = []
    #params_true = np.array([m1/MSOLAR, m2/MSOLAR])
    params_true = np.array([getChirpMass(m1, m2)/MSOLAR, getTotalMass(m1, m2)/MSOLAR])
    for i in range(len(beta)):
        #starting guess for masses
        m1_guess = 0.6*MSOLAR
        m2_guess = 0.6*MSOLAR
        fdot = beta[i]/(t_obs**2)
        chirp_guess = (fdot * 5 / (96 * (np.pi**(8/3)) * (freq0**(11/3))))**(3./5)
        total_guess = chirp_guess / (0.24)**(3./5)
        #final_guess = getRootFinder_tides_componentMass(freq0, beta[i], delta[i], t_obs, params_true[0], params_true[1], m1_guess, m2_guess)
        final_guess = getRootFinder_tides_chirpTotalMass(freq0, beta[i], delta[i], t_obs, params_true[0], params_true[1], chirp_guess, total_guess)
        #filter the masses
        if final_guess.success == True:
            if np.isnan(final_guess.x[0]) == False and np.isnan(final_guess.x[1]) == False:
                chirpMass.append(final_guess.x[0])
                totalMass.append(final_guess.x[1])
    #now covert to m1, m2
    makeHistogramPlots(chirpMass, "Mc")
    makeHistogramPlots(totalMass, "Mt")
    convertedM1 = []
    convertedM2 = []
    for i in range(len(chirpMass)):
        masses = convertChirpTotal_to_M1M2(chirpMass[i], totalMass[i])
        if masses[0] < 1.4 and masses[1] < 1.4:
            convertedM1.append(masses[0])
            convertedM2.append(masses[1])

    # print("Real Mc:", params_true[0])
    # print("final mc:", np.mean(chirpMass))
    # print("real Mt:", params_true[1])
    # print("final mt:", np.mean(totalMass))

    # print("M1:", np.mean(convertedM1))
    # print("M2:", np.mean(convertedM2))

    return convertedM1, convertedM2

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