
from matplotlib import pyplot as plt
from numpy import number
import numpy as np
import matplotlib.ticker
from scipy.optimize import curve_fit
from globals import *

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
    mass_range = []
    R_range = []
    for i in range(len(masses)):
        if masses[i] >= 0.1:
            mass_range.append(masses[i])
            R_range.append(r[i])
    rm = 10**9 * (np.asarray(mass_range) / 0.6 )**(-1/3)

    #get the polynomial fit for the MvR relation 
    # poly = np.polyfit(np.log(masses), np.log(r), 1)
    # poly_simple = np.polyfit(np.log(masses), np.log(rm), 1)
    # print("Numerical result:",poly)
    # print("Analytic result:",poly_simple)   

    fig1, ax1 = plt.subplots()
    plt.rcParams.update({'font.size': 19})
    ax1.loglog(R_range, mass_range, label="Numerical", linewidth=2)
    ax1.loglog(rm, mass_range, label="Kuns Scaling Relation", linewidth=2)
    
    ax1.set_yticks([0.1,0.5, 1])
    ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax1.set_xticks([10**9,2*10**9])
    ax1.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())    
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
        if masses[i] >= 0.1:
            mass_range.append(masses[i])
            I0_range.append(I0[i])

    Im = 3.1e50 * (np.asarray(mass_range) /0.6)**(1/3)
    #get the polynomial fit for the MvR relation 
    # poly = np.polyfit(np.log(masses), np.log(I0), 1)
    # poly_cut = np.polyfit(np.log(mass_range), np.log(I0_range),1)
    # poly_simple = np.polyfit(np.log(masses), np.log(Im), 1)
    # print("Numerical result:",poly)
    # print("Numerical Fixed result:", poly_cut)
    # print("Analytic result:",poly_simple)   

    fig1, ax1 = plt.subplots()
    plt.rcParams.update({'font.size': 19})
    ax1.loglog(mass_range, I0_range, label="Numerical", linewidth=2)
    #ax1.loglog(mass_range, I0_range, label="Numerical", linewidth=2)
    ax1.loglog(mass_range, Im, label="Kuns Scaling Relation", linewidth=2)

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

def plot3DMassRadiusVelocity(t, totalMasses, omega):
    ax = plt.figure().add_subplot(projection='3d')

    totalMasses.sort()
    x = np.array(omega)
    y = np.array(t)
    z = np.array(totalMasses)

    X, Y = np.meshgrid(x, y)

    # Plot the 3D surface
    ax.plot_surface(X, Y, z, edgecolor='royalblue')
    ax.set_xlabel('Rotational Velocity', labelpad=20)
    ax.set_ylabel('Radius', labelpad=20)
    ax.set_zlabel('Mass', labelpad=20)
    ax.set_title("Mass-Radius-Spin Relation")
    plt.show()

def plot3DMassInertiaVelocity(totalMasses, totalInertia, omega):
    ax = plt.figure().add_subplot(projection='3d')

    totalMasses.sort()
    totalInertia.sort()
    x = np.array(omega)
    y = np.array(totalMasses)
    z = np.array(totalInertia)

    # Plot the 3D surface
    ax.plot_surface(x, y, z, edgecolor='royalblue')
    ax.set_xlabel('Rotational Velocity - x', labelpad=20)
    ax.set_ylabel('Mass - y', labelpad=20)
    ax.set_zlabel('Inertia - z', labelpad=20)
    ax.set_title("Mass-Inertia-Spin Relation")
    plt.show()

def plot2DMassInertiaVelocity(totalMasses, totalInertia):
    #need to be 2D arrays. We want the ROWS of the inputs
    for i in range(len(totalInertia)):
        totalMasses[i:].sort()
        totalInertia[i:].sort()

    #we're taking second to last so it's not 0
    MassN = [mass[-2] for mass in totalMasses]
    inertiaN = [inertia[-2] for inertia in totalInertia]

    ref_masses = []
    ref_inertias = []

    for i in range(len(MassN)):
        if MassN[i] >= 0.20:
            ref_masses.append(MassN[i])
            ref_inertias.append(inertiaN[i]/10**48)

    #fitting the scaling factors for range of spins from above
    alpha_alt = [1,(19.85*omg)**2, (13.2*omg)**2, (9.9*omg)**2, (7.9*omg)**2, (6.6*omg)**2, 
                  (5.65*omg)**2, (4.95*omg)**2, (4.4*omg)**2, (3.95*omg)**2]
    omega = np.linspace(0.1*omg, omg, 10)
    
    #polynomial fitting function
    def f1(x, a0, a1, a2):
        return a0 + a1/x + a2/np.power(x,2)
   
    popt1, pcov1 = curve_fit(f1, ref_masses, ref_inertias)
    popt2, pcov2 = curve_fit(f1, omega, alpha_alt)

    IvsM_omega_fit = f1(ref_masses, popt1[0], popt1[1], popt1[2])

    # print( "Parameters from I(2) vs M fit:")
    # print( "a =", popt1[0], "+/-", pcov1[0,0]**0.5)
    # print ("b =", popt1[1], "+/-", pcov1[1,1]**0.5)
    # print("c =", popt1[2], "+/-", pcov1[2,2]**0.5)
    # print( "Parameters from alpha_alt fit:")
    # print( "a =", popt2[0], "+/-", pcov2[0,0]**0.5)
    # print ("b =", popt2[1], "+/-", pcov2[1,1]**0.5)
    # print("c =", popt2[2], "+/-", pcov2[2,2]**0.5)

    plt.figure()
    for i in range(len(omega)):  
        alpha_polynomial_fit = f1(omega[i], popt2[0], popt2[1], popt2[2])
        total_fit = IvsM_omega_fit / alpha_polynomial_fit   
        scale = round(omega[i] / omg, 1) 
        plt.plot(ref_masses, total_fit, label = "total fit, $\Omega$ = %s$\Omega_{max}$"%scale)

    plt.title("Total Fit for Moment of Inertia and Mass as Function of Spin")
    plt.xlabel("Total Mass")
    plt.ylabel("Moment of Inertia, $I^{(2)}$")
    plt.yscale("log")
    plt.legend()
    plt.show()