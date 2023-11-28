import numpy as np
from const import *

def getChirpMass(mass1, mass2):
    return(((mass1*mass2)**(3/5)) / ((mass1 + mass2)**(1/5)))

def getTotalMass(mass1, mass2):
    return(mass1 + mass2)

def convertChirpTotal_to_M1M2(chirpmass, totalmass):
    m1 = (1/2) * (totalmass + np.sqrt(totalmass**2 - (4*chirpmass**(5/3)*totalmass**(1/3))))
    m2 = (1/2) * (totalmass - np.sqrt(totalmass**(1/3) * (totalmass**(5/3) - 4*chirpmass**(5/3))))
    q = m2 / m1
    return [m1,m2,q]

def get_comp_mass_Mc_Mt(Mc, Mt):
    eta = (Mc/Mt)**(5./3.)
    qq = ((1-2.*eta) - np.sqrt((1-2.*eta)**2. - 4.*eta**2.))/(2.*eta)
    M1=Mt/(1.+qq)
    M2=qq*Mt/(1.+qq)
    q = M2 / M1
    return M1, M2, q 

def calculateAmplitude(SNR, t_obs):
    #calculates approx amplitude from SNR^2 = int( h(t)^2 dt) - time avg over 0 to t_obs
    return(SNR / np.sqrt(t_obs))

def calculateAmplitude_phys(dl, chirpMass, freq0):
    #calculates approx amplitude from SNR^2 = int( h(t)^2 dt) - time avg over 0 to t_obs
    return(np.pi**2/3 * chirpMass**(5/3) * freq0**(2/3) / dl)

def derivative(fx, params, h, i = 0):
    np.isscalar(params)
    isParamsArray = not np.isscalar(params)
    isHArray = not np.isscalar(h)

    if isParamsArray != isHArray:
        print("ERROR: Params and h must be the same size")
        return "ERROR"
    elif isParamsArray != True:
        params = [params]
        h = [h]
    elif len(params) != len(h):
        print("ERROR: Params and h must be the same size")
        return "ERROR"

    return (fx(getParamsWithStep(params, i, h[i])) - fx(getParamsWithStep(params, i, h[i], False))) / (2*h[i])

def getParamsWithStep(params, target, step, stepUp = True):
    newParams = np.copy(params)

    if (stepUp):
        newParams[target] = params[target] + step
    else:
        newParams[target] = params[target] - step

    return newParams

import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq

def dataFFT(freq_file, time_file, duration, samplerate):
    #import signal file
    signal = []
    times = []
    f = open(time_file,'r') 
    freqfile = open(freq_file, 'r')
    for row in f: 
        row = row.split('\n') 
        times.append(float(row[0]))  
    for row in freqfile: 
        row = row.split('\n') 
        signal.append(float(row[0]))

    signal_cut = []
    times_cut = []
    for sig in range(len(signal)):
        if times[sig] <= duration:
            signal_cut.append(signal[sig])
            times_cut.append(times[sig])

    N = len(signal_cut)
    h_f = fft(signal_cut)
    f_t = fftfreq(N, 1/samplerate)

    import PhenomA as pa
    import LISA as li
    # create LISA object
    lisa = li.LISA() 

    # Plot LISA's sensitivity curve
    f  = np.logspace(np.log10(1.0e-5), np.log10(1.0e0), N)
    Sn = lisa.Sn(f)

    plt.figure()
    plt.plot(times_cut, signal_cut)
    plt.show()
    plt.figure()
    plt.loglog(f, np.abs(h_f))
    plt.loglog(f, np.sqrt(f*Sn))
    plt.xlabel("Freq (Hz)")
    plt.ylabel("h(f)")
    plt.xlim(1.0e-5, 1.0e0)
    plt.ylim(3.0e-22, 1.0e-15)
    plt.show()

    return
    

