import numpy as np
from const import *
from graphing import * 
from postprocessing import betaDeltaM1M2Converter

def create_prior_model(params, sigma_prior_lim, sigmas):
    """get the limits for the prior model, given fiducial starting parameters,
    an estimate of the standard deviations (only matters for amp0,freq0,freqD currently),
    and how many standard deviations to include (subject to other constraints)
    inputs:
        params: 1D float array, fiducial parameters"""
    #define the parameter indicies
    idx_beta = 0
    idx_delta = 1

    n_par = np.size(params)
    low_lims = np.zeros(n_par)
    high_lims = np.zeros(n_par)

    low_lims[idx_beta] = params[idx_beta]-sigma_prior_lim*sigmas[idx_beta]
    high_lims[idx_beta] = params[idx_beta]+sigma_prior_lim*sigmas[idx_beta]
    
    low_lims[idx_delta] = params[idx_delta]-sigma_prior_lim*sigmas[idx_delta]
    high_lims[idx_delta] = params[idx_delta]+sigma_prior_lim*sigmas[idx_delta]

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
            fx = lambda params: Frequency_Tides_MassRatio(freq0, params, t_obs)[i]
            dF = lambda x, params: derivative(fx, params, step_size, x)
            jacobian[i, j] = dF(j, p)
    #invert this matrix, need d(m1m2)/d(b,d)
    inverseJac = np.linalg.inv(jacobian)
    #take determinant of jacobian
    detjac = np.abs(np.linalg.det(inverseJac))
    return detjac

def resampling(params, freq0, t_obs, nsteps, sigmas):
    #establish prior:
    sigma_prior_lim=2
    #high_lims, low_lims = create_prior_model(params, sigma_prior_lim, sigmas)
    masses = [0.7*MSOLAR, 0.6*MSOLAR]

    #determine m1, m2 for high and low limits
    #beta_low, delta_low = Frequency_Tides_Masses(freq0, [0.4*MSOLAR, 0.4*MSOLAR], t_obs)
    #beta_up, delta_up = Frequency_Tides_Masses(freq0, [1.4*MSOLAR, 1.4*MSOLAR], t_obs)
    high_lims = [4000., 3.2]
    low_lims = [1000., 0.8]
    #generate draws and convert to beta, delta
    beta_prior = []
    delta_prior = []
    
    i=0
    while i < nsteps:
        #draw from prior - returns array of the parameter draws
        prior = prior_draw(low_lims, high_lims, params)
        beta_prior.append(prior[0])
        delta_prior.append(prior[1])
        i += 1
    print(i)
    #now run these through to get m1, m2
    mass1_prior = []
    mass2_prior = []
    massratio_prior = []

    jacobian = []
    mass1, mass2, q = betaDeltaM1M2Converter(beta_prior, delta_prior, freq0, t_obs, masses[0], masses[1])
    for i in range(len(mass1)):
        if mass1[i] >= 0. and mass1[i] <= 1.4:
            if q[i] >= 0.:
                mass1_prior.append(mass1[i])
                mass2_prior.append(mass2[i])
                massratio_prior.append(q[i])
                j = get_Jacobian([mass1[i]*MSOLAR, q[i]], freq0, t_obs)
                jacobian.append(j)
    print(len(mass1_prior))
    makeCornerPlot([beta_prior, delta_prior], [np.average(beta_prior), np.average(delta_prior)], labels=["beta", "delta"])

    masses_Msun = [0.7, 0.857143]
    makeCornerPlot([mass1_prior, massratio_prior], masses_Msun, weight=jacobian, labels=[r"$M_{1}$", "q"])

    return 

def drawSamples_M1M2(nsteps, low_lims, high_lims, params, freq0):
    m1_prior = []
    m2_prior = []
    i=0
    while i < nsteps:
        prior = prior_draw(high_lims, low_lims, params)
        m1_prior.append(prior[0])
        m2_prior.append(prior[1])
        i += 1
    #get beta, delta
    beta_prior = []
    delta_prior = []
    for j in range(len(m1_prior)):
        beta, delta = Frequency_Tides_Masses(freq0, [m1_prior[j]*MSOLAR, m2_prior[j]*MSOLAR], t_obs=4*SECSYEAR)
        beta_prior.append(beta)
        delta_prior.append(delta)
    makeCornerPlot([beta_prior, delta_prior, m1_prior, m2_prior], [3600, 1.38, 0.485, 0.485], labels=["beta", "delta", "M1", "M2"])
    return



