import numpy as np
from const import *
from graphing import * 

def create_prior_model(params):
    """get the limits for the prior model, given fiducial starting parameters,
    an estimate of the standard deviations (only matters for amp0,freq0,freqD currently),
    and how many standard deviations to include (subject to other constraints)
    inputs:
        params: 1D float array, fiducial parameters"""
    #define the parameter indicies
    idx_mass1 = 0
    idx_mass2 = 1

    n_par = np.size(params)
    low_lims = np.zeros(n_par)
    high_lims = np.zeros(n_par)

    low_lims[idx_mass1] = 0.
    high_lims[idx_mass1] = 1.4*MSOLAR
    
    low_lims[idx_mass2] = 0.
    high_lims[idx_mass2] = 1.4*MSOLAR

    return high_lims, low_lims

def prior_draw(low_lims, high_lims, params):
    """do a prior draw"""
    n_par = np.size(params)
    draw = np.zeros(n_par)
    for itrp in range(n_par):
        draw[itrp] = np.random.uniform(low_lims[itrp], high_lims[itrp])
    return draw

def get_Jacobian(p, freq0, t_obs):
    step_size = [1.e-10, 1.e-10]
    jacobian = np.zeros((np.size(p), np.size(p)))
    for i in range(len(p)):
        for j in range(len(p)):
            fx = lambda params: Frequency_Tides_Masses(freq0, params, t_obs)[i]
            dF = lambda x, params: derivative(fx, params, step_size, x)
            jacobian[i, j] = dF(j, p)
    #take determinant of jacobian
    detjac = np.abs(np.linalg.det(jacobian))
    return detjac

def resampling(params, freq0, t_obs, nsteps):
    #establish prior:
    high_lims, low_lims = create_prior_model(params)
    #establish jacobian - returns abs value of det(J)
    jacobian = get_Jacobian(params, freq0, t_obs)
    #generate draws and convert to beta, delta
    freq_prior = []
    m1_prior = []
    m2_prior = []
    
    i=0
    while i < nsteps:
        #draw from prior - returns array of the parameter draws
        prior = prior_draw(low_lims, high_lims, params)
        m1_prior.append(prior[0])
        m2_prior.append(prior[1])
        weight = prior[0] * prior[1] * jacobian
        freq_prior.append(weight)
        i += 1
    print(i)
    makeHistogramPlots(freq_prior, "freq prior")
    makeHistogramPlots(m1_prior, "M1")
    makeHistogramPlots(m2_prior, "M2")
    return freq_prior, m1_prior, m2_prior

