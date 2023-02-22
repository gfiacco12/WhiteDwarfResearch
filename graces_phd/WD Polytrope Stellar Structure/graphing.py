
from matplotlib import pyplot as plt
from numpy import number

# plotting central density vs final mass


def plotCentralDensity(masses: "list[number]", rhos: "list[number]"):
    plt.figure()
    plt.plot(rhos, masses)
    plt.title("Mass v Central Density for Polytropic EOS")
    plt.xlabel("Central Density $\\rho_{0}$")
    plt.ylabel("Final Mass ($M_{\\odot}$)")
    plt.show()

def plotMassRadiusRelation(masses: "list[number]", r: "list[number]"):
    plt.figure()
    plt.plot(r, masses)
    plt.title("Mass v Radius Relationship for Polytropic WD")
    plt.xlabel("Radius (cm)")
    plt.ylabel("Final Mass ($M_{\\odot}$)")
    plt.show()

def plotStellarStructure(t, densities: "list[number]", totalMasses: "list[number]", totalPressures: "list[number]"):
    # plotting total quantities vs R - including second order corrections
    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    ax1.plot(t, totalMasses)
    ax1.set_title("M vs r")
    ax1.set_xlabel("Radius (cm)")
    ax1.set_ylabel("Total Mass ($M_{\\odot})$")
    ax2.plot(t, densities)
    ax2.set_title("$\\rho$ vs r")
    ax2.set_ylabel('Density (g/cm^3)')
    ax2.set_xlabel('Radius (cm)')
    ax3.plot(t, totalPressures)
    ax3.set_title("p vs r")
    ax3.set_ylabel('Pressure')
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