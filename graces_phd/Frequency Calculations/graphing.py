import numpy as np
import matplotlib.pyplot as plt
from const import *
from FrequencyCalculations import *
from FisherCalculations import *
import corner

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
    print(delta_tide_list)
    plt.figure()
    plt.semilogy(freq_range, delta_tide_list, label='Tides')
    plt.semilogy(freq_range, delta_1PN_list, label='1PN')
    plt.semilogy(freq_range, fisher_errors, label="$\Delta\delta$")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("$\delta$")
    plt.title("Relationship between $\delta$ and Frequency")
    plt.legend()
    plt.show()

def makeHistogramPlots(data_file, title, lowlim, uplim):

    f = open(data_file, 'r')
    massData = []
    for row in f:
        #elements = row.split(' ')
        elements = row.split(' ')
        elements = list(map(lambda e : float(e), elements))
        massData += elements

    plt.figure()
    plt.hist(massData, 20)
    plt.xlim(lowlim, uplim)
    plt.title(title)
    plt.show()

def makeCornerPlot(chirpmass, totalmass, real_params):
    
    f = open(chirpmass, 'r')
    g = open(totalmass, 'r')
    chirpMass = []
    totalMass = []
    for row in f:
        #elements = row.split(' ')
        elements = row.split(' ')
        elements = list(map(lambda e : float(e), elements))
        chirpMass += elements
    for row in g:
        #elements = row.split(' ')
        elements = row.split(' ')
        elements = list(map(lambda e : float(e), elements))
        totalMass += elements
    samples = np.array([chirpMass, totalMass])

    params_true = real_params
    labels = [r"$\mathcal{M}_{c}$", r"$M_{T}$"]
    # create the corner plot figure
    fig = plt.figure(figsize=(10, 7.5))
    figure = corner.corner(samples.T, fig=fig, bins=25, hist_kwargs={"density": True}, show_titles=True, title_fmt=None,
                           title_kwargs={"fontsize": 12}, labels=labels, max_n_ticks=3, label_kwargs={"fontsize": 12}, labelpad=0.15,
                           smooth=0.25, title_quantiles=[0.16,0.5,0.84], quantiles= [0.16,0.5,0.84], levels=[0.682, 0.954], truths=params_true)
    # overplot the true parameters
    corner.overplot_points(figure, params_true[None], marker="s", color='tab:blue', markersize=4)
    corner.overplot_lines(figure, params_true, color='tab:blue')
    # adjust the figure to fit the box better
    fig.subplots_adjust(wspace=0., hspace=0., left=0.05, top=0.95, right=0.99, bottom=0.05)
    for ax in figure.get_axes():
        ax.tick_params(which='both', direction='in', bottom=True, top=True, left=True, right=True, labelsize=6)
    plt.show()