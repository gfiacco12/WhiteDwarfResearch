import numpy as np
import matplotlib.pyplot as plt
from const import *
from FrequencyCalculations import *
from FisherCalculations import *

def makeDeltaPlot(freq0, mass1, mass2, tobs, amp):

    #make range of frequencies from 10-20mHz
    freq_range = np.linspace(10.e-3, freq0, 10)
    fisher_errors = np.linspace(2.3*10**-2, 2.3*10**-2, 10)
    delta_tide_list = []
    delta_1PN_list = []
 
    #calculate delta for both cases for each freq
    for i in freq_range:
        freq_1PN = Frequency_1PN(i, mass1, mass2, tobs)
        freq_tide = Frequency_Tides_Masses(i, mass1, mass2, tobs)

        delta_1PN = freq_1PN[2] * tobs**3
        delta_tide = freq_tide[2] * tobs**3

        delta_1PN_list.append(delta_1PN)
        delta_tide_list.append(delta_tide)
    plt.figure()
    plt.semilogy(freq_range, delta_tide_list, label='Tides')
    plt.semilogy(freq_range, delta_1PN_list, label='1PN')
    plt.semilogy(freq_range, fisher_errors, label="$\Delta\delta$")
    plt.xlabel("Frequency (mHz)")
    plt.ylabel("$\delta$")
    plt.title("Relationship between $\delta$ and Frequency")
    plt.legend()
    plt.show()
