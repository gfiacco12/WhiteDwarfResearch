
from matplotlib import pyplot as plt
from numpy import number
import numpy as np
import matplotlib.ticker

# plotting central density vs final mass


def plotCentralDensity(masses: "list[number]", rhos: "list[number]"):
    masses.sort()
    rhos.sort()
    plt.figure()
    plt.plot(rhos, masses)
    plt.title("Mass v Central Density for Polytropic EOS")
    plt.xlabel("Central Density $\\rho_{0}$ ($g/cm^{3}$)")
    plt.ylabel("Final Mass ($M_{\\odot}$)")
    plt.show()

def plotMassRadiusRelation(masses: "list[number]", r: "list[number]"):
    #plot the scaling relation from Kuns 2020
    masses.sort()
    rm = 10**9 * (np.asarray(masses) / 0.6 )**(-1/3)

    #get the polynomial fit for the MvR relation 
    poly = np.polyfit(np.log(masses), np.log(r), 1)
    poly_simple = np.polyfit(np.log(masses), np.log(rm), 1)
    print("Numerical result:",poly)
    print("Analytic result:",poly_simple)   

    fig1, ax1 = plt.subplots()
    plt.rcParams.update({'font.size': 19})
    ax1.loglog(r, masses, label="Numerical", linewidth=2)
    ax1.loglog(rm, masses, label="Kuns Scaling Relation", linewidth=2)
    
    ax1.set_yticks([0.1,0.5, 1])
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())    
    ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')

    plt.title("Mass v Radius Relationship for Polytropic WD")
    plt.xlabel("Radius (cm)", fontsize=16)
    plt.ylabel("Final Mass ($M_{\\odot}$)",fontsize=16)
    plt.legend(fontsize=15,loc='upper right')
    plt.show()

def plotMassInertiaRelation(masses: "list[number]", I0: "list[number]"):
    masses.sort()
    I0.sort() 
    mass_range = []
    I0_range = []
    for i in range(len(masses)):
        if masses[i] >= 0.29:
            mass_range.append(masses[i])
            I0_range.append(I0[i])

    Im = 3.1e50 * (np.asarray(masses) /0.6)**(1/3)
    #get the polynomial fit for the MvR relation 
    poly = np.polyfit(np.log(masses), np.log(I0), 1)
    poly_cut = np.polyfit(np.log(mass_range), np.log(I0_range),1)
    poly_simple = np.polyfit(np.log(masses), np.log(Im), 1)
    print("Numerical result:",poly)
    print("Numerical Fixed result:", poly_cut)
    print("Analytic result:",poly_simple)   

    fig1, ax1 = plt.subplots()
    plt.rcParams.update({'font.size': 19})
    ax1.loglog(masses, I0, label="Numerical", linewidth=2)
    #ax1.loglog(mass_range, I0_range, label="Numerical", linewidth=2)
    ax1.loglog(masses, Im, label="Kuns Scaling Relation", linewidth=2)

    ax1.set_xticks([0.1,0.5, 1])
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())    
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.set_aspect(1.0/ax1.get_data_ratio(), adjustable='box')

    plt.title("Mass v $I^{(0)}$ Relationship for Polytropic WD")
    plt.xlabel("Mass ($M_{\\odot}$)",fontsize=16)
    plt.ylabel("$I^{(0)}$ ($g cm^{2}$)",fontsize=16)
    plt.legend()
    plt.show()

def plotStellarStructure(t, densities: "list[number]", totalMasses: "list[number]", totalPressures: "list[number]"):
    # plotting total quantities vs R - including second order corrections
    plt.rcParams.update({'font.size': 20})
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    fig.tight_layout(pad=0.75)
    ax1.plot(t, totalMasses)
    ax1.set_title("M vs r")
    ax1.set_xlabel("Radius (cm)")
    ax1.set_ylabel("Total Mass ($M_{\\odot})$")
    ax2.plot(t, densities)
    ax2.set_title("$\\rho$ vs r")
    ax2.set_ylabel('Density ($g cm^{-3}$)')
    ax2.set_xlabel('Radius (cm)')
    ax3.plot(t, totalPressures)
    ax3.set_title("p vs r")
    ax3.set_ylabel('Pressure ($g cm^{-1} s^{-2}$)')
    ax3.set_xlabel('Radius (cm)')
    plt.show()

def plotMomentofInertia(t, I0: "list[number]", I2: "list[number]", totalI: "list[number]"):
    #plotting the different moments of inertia
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(t, I0)
    ax1.set_title("$I^{0}$ vs r")
    ax1.set_xlabel("Radius (cm)")
    ax1.set_ylabel("Oth Order I (g*cm^2)")
    ax2.plot(t, I2)
    ax2.set_title("$I^{2}$ vs r")
    ax2.set_ylabel('Rotational Correction to I (g*cm^2)')
    ax2.set_xlabel('Radius (cm)')
    ax3.plot(t, totalI)
    ax3.set_title("Total Moment of Inertia vs r")
    ax3.set_ylabel('Total I (g*cm^2)')
    ax3.set_xlabel('Radius (cm)')
    plt.show()