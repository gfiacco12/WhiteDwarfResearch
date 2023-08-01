
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
    
    Mass0 = [mass[0] for mass in totalMasses]
    Mass1 = [mass[1] for mass in totalMasses]
    Mass2 = [mass[2] for mass in totalMasses]
    Mass3 = [mass[3] for mass in totalMasses]
    Mass4 = [mass[4] for mass in totalMasses]
    Mass5 = [mass[5] for mass in totalMasses]
    Mass6 = [mass[6] for mass in totalMasses]
    Mass7 = [mass[7] for mass in totalMasses]
    Mass8 = [mass[8] for mass in totalMasses]
    Mass9 = [mass[9] for mass in totalMasses]
    Mass10 = [mass[10] for mass in totalMasses]

    inertia0 = [inertia[0] for inertia in totalInertia]
    inertia1 = [inertia[1] for inertia in totalInertia]
    inertia2 = [inertia[2] for inertia in totalInertia]
    inertia3 = [inertia[3] for inertia in totalInertia]
    inertia4 = [inertia[4] for inertia in totalInertia]
    inertia5 = [inertia[5] for inertia in totalInertia]
    inertia6 = [inertia[6] for inertia in totalInertia]
    inertia7 = [inertia[7] for inertia in totalInertia]
    inertia8 = [inertia[8] for inertia in totalInertia]
    inertia9 = [inertia[9] for inertia in totalInertia]
    inertia10 = [inertia[10] for inertia in totalInertia]

    Mass_max = []
    Mass_1 = []
    Mass_2 = []
    Mass_3 = []
    Mass_4 = []
    Mass_5 = []
    Mass_6 = []
    Mass_7 = []
    Mass_8 = []
    Mass_9 = []
    Mass_10 = []
    inertia_max = []
    inertia_1 = []
    inertia_2 = []
    inertia_3 = []
    inertia_4 = []
    inertia_5 = []
    inertia_6 = []
    inertia_7 = []
    inertia_8 = []
    inertia_9 = []
    inertia_10 = []

    for i in range(len(Mass0)):
        if Mass0[i] >= 0.205:
            Mass_max.append(Mass0[i])
            inertia_max.append(inertia0[i]/10**48)
        if Mass1[i] >= 0.205:
            Mass_1.append(Mass1[i])
            inertia_1.append(inertia1[i]/10**48)
        if Mass2[i] >= 0.205:
            Mass_2.append(Mass2[i])
            inertia_2.append(inertia2[i]/10**48)
        if Mass3[i] >= 0.205:
            Mass_3.append(Mass3[i])
            inertia_3.append(inertia3[i]/10**48)
        if Mass4[i] >= 0.202:
            Mass_4.append(Mass4[i])
            inertia_4.append(inertia4[i]/10**48)
        if Mass5[i] >= 0.202:
            Mass_5.append(Mass5[i])
            inertia_5.append(inertia5[i]/10**48)
        if Mass6[i] >= 0.201:
            Mass_6.append(Mass6[i])
            inertia_6.append(inertia6[i]/10**48)
        if Mass7[i] >= 0.201:
            Mass_7.append(Mass7[i])
            inertia_7.append(inertia7[i]/10**48)
        if Mass8[i] >= 0.201:
            Mass_8.append(Mass8[i])
            inertia_8.append(inertia8[i]/10**48)
        if Mass9[i] >= 0.2:
            Mass_9.append(Mass9[i])
            inertia_9.append(inertia9[i]/10**48)
        if Mass10[i] >= 0.2:
            Mass_10.append(Mass10[i])
            inertia_10.append(inertia10[i]) #this one is 0, does not show up

    #fitting the scaling factors for range of spins from above
    alpha_alt = [1,(19.85*omg)**2, (13.2*omg)**2, (9.9*omg)**2, (7.9*omg)**2, (6.6*omg)**2, 
                  (5.65*omg)**2, (4.95*omg)**2, (4.4*omg)**2, (3.95*omg)**2]
    omega = np.linspace(0.1*omg, omg, 10)
    
    #polynomial fitting function
    def f1(x, a0, a1, a2):
        return a0 + a1/x + a2/np.power(x,2)

    #fit the product of the polys
    omega_range = np.linspace(min(omega), max(omega), len(Mass_9))
    mass_range = np.linspace(min(Mass_9), max(Mass_9), len(Mass_9))
    mega_mass_list = np.concatenate([Mass_9, Mass_8, Mass_7, Mass_6, Mass_5, Mass_4, Mass_3, Mass_2, Mass_1, Mass_max])
    mega_omega_list = np.concatenate([omega[0]*np.ones(len(Mass_9)), omega[1]*np.ones(len(Mass_8)), omega[2]*np.ones(len(Mass_7)), omega[3]*np.ones(len(Mass_6)),
                                      omega[4]*np.ones(len(Mass_5)), omega[5]*np.ones(len(Mass_4)), omega[6]*np.ones(len(Mass_3)),
                                      omega[7]*np.ones(len(Mass_2)), omega[8]*np.ones(len(Mass_1)), omega[9]*np.ones(len(Mass_max))])
    mega_inertia_list = np.concatenate([inertia_9, inertia_8, inertia_7, inertia_6, inertia_5, inertia_4, inertia_3, inertia_2, inertia_1, inertia_max])
   
    popt1, pcov1 = curve_fit(f1, Mass_9, inertia_9)
    popt2, pcov2 = curve_fit(f1, omega, alpha_alt)
    print( "Parameters from I(2) vs M fit:")
    print( "a =", popt1[0], "+/-", pcov1[0,0]**0.5)
    print ("b =", popt1[1], "+/-", pcov1[1,1]**0.5)
    print("c =", popt1[2], "+/-", pcov1[2,2]**0.5)
    print( "Parameters from alpha_alt fit:")
    print( "a =", popt2[0], "+/-", pcov2[0,0]**0.5)
    print ("b =", popt2[1], "+/-", pcov2[1,1]**0.5)
    print("c =", popt2[2], "+/-", pcov2[2,2]**0.5)

    #create a total fit - normalize the alpha fit to 1 at omega_max
    alpha_poly_ref = f1(0.1*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin2 = f1(0.2*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin3 = f1(0.3*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin4 = f1(0.4*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin5 = f1(0.5*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin6 = f1(0.6*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin7 = f1(0.7*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin8 = f1(0.8*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin9 = f1(0.9*omg, popt2[0], popt2[1], popt2[2])
    alpha_poly_spin_max = f1(omg, popt2[0], popt2[1], popt2[2])
    IvsM_omega_fit = f1(Mass_9, popt1[0], popt1[1], popt1[2])
    ref_total_fit =  (IvsM_omega_fit)*(alpha_poly_ref)
    spin2_total_fit = (IvsM_omega_fit)/(alpha_poly_spin2)
    spin3_total_fit = (IvsM_omega_fit)/(alpha_poly_spin3)
    spin4_total_fit = (IvsM_omega_fit)/(alpha_poly_spin4)
    spin5_total_fit = (IvsM_omega_fit)/(alpha_poly_spin5)
    spin6_total_fit = (IvsM_omega_fit)/(alpha_poly_spin6)
    spin7_total_fit = (IvsM_omega_fit)/(alpha_poly_spin7)
    spin8_total_fit = (IvsM_omega_fit)/(alpha_poly_spin8)
    spin9_total_fit = (IvsM_omega_fit)/(alpha_poly_spin9)
    spin_max_total_fit = (IvsM_omega_fit)/(alpha_poly_spin_max)


    #plotting I(2) vs M
    plt.figure()
    plt.plot(Mass_max, inertia_max, label="$\Omega_{max}$")
    plt.plot(Mass_1, inertia_1, label="$\Omega$ = 0.9$\Omega_{max}$")
    plt.plot(Mass_2, inertia_2, label="$\Omega$ = 0.8$\Omega_{max}$")
    plt.plot(Mass_3, inertia_3, label="$\Omega$ = 0.7$\Omega_{max}$")
    plt.plot(Mass_4, inertia_4, label="$\Omega$ = 0.6$\Omega_{max}$")
    plt.plot(Mass_5, inertia_5, label="$\Omega$ = 0.5$\Omega_{max}$")
    plt.plot(Mass_6, inertia_6, label="$\Omega$ = 0.4$\Omega_{max}$")
    plt.plot(Mass_7, inertia_7, label="$\Omega$ = 0.3$\Omega_{max}$")
    plt.plot(Mass_8, inertia_8, label="$\Omega$ = 0.2$\Omega_{max}$")
    plt.plot(Mass_9, inertia_9, label="$\Omega$ = 0.1$\Omega_{max}$")
    #plt.plot(Mass_10, inertia_10, label="$\Omega$ = 0")
    plt.yscale("log")
    plt.xlabel("Total Mass")
    plt.ylabel("2nd Order Moment of Inertia ($x10^{48}$)")
    plt.title("Mass vs Moment of Inertia for Range of Spin")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(Mass_9, ref_total_fit, label="total fit, $\Omega$ = 0.1$\Omega_{max}$")
    plt.plot(Mass_8, spin2_total_fit, label="total fit, $\Omega$ = 0.2$\Omega_{max}$")
    plt.plot(Mass_7, spin3_total_fit, label="total fit, $\Omega$ = 0.3$\Omega_{max}$")
    plt.plot(Mass_6, spin4_total_fit, label="total fit, $\Omega$ = 0.4$\Omega_{max}$")
    plt.plot(Mass_5, spin5_total_fit, label="total fit, $\Omega$ = 0.5$\Omega_{max}$")
    plt.plot(Mass_4, spin6_total_fit, label="total fit, $\Omega$ = 0.6$\Omega_{max}$")
    plt.plot(Mass_3, spin7_total_fit, label="total fit, $\Omega$ = 0.7$\Omega_{max}$")
    plt.plot(Mass_2, spin8_total_fit, label="total fit, $\Omega$ = 0.8$\Omega_{max}$")
    plt.plot(Mass_1, spin9_total_fit, label="total fit, $\Omega$ = 0.9$\Omega_{max}$")
    plt.plot(Mass_max, spin_max_total_fit, label="total fit, $\Omega$ = $\Omega_{max}$")
    plt.title("Total Fit for Moment of Inertia and Mass as Function of Spin")
    plt.xlabel("Total Mass")
    plt.ylabel("Moment of Inertia, $I^{(2)}$")
    plt.yscale("log")
    plt.legend()
    plt.show()