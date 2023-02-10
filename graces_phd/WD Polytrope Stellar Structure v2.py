# MHD Homework 2b - Integrating Hydrostatic Equations - Grace Fiacco
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import integrate
from scipy import optimize

# start by defining global variables needed for all the equations - in cgs ick
# Constants needed: G, mp, mu, kb, sig_sb,hbar, electron mass
G = 6.67e-8
MP = 1.6726e-24  # g
MU = 0.8
KB = 1.38e-16
SIG = 5.67e-5
H = 6.626e-27
ME = 9.1094e-28
K = (H ** 2 / (5*ME)) * (3/(8*np.pi))**(2/3)
omg = (G*1.193e33/((5.5e9)**3))**(1/2)  # angular velocity - keplerian


def hydro_ODE(r, y):
    rho = y[0]
    M = y[1]
    p = y[2]
    # l=0
    p2 = y[3]
    M2 = y[4]
    # l=2, particular
    X = y[5]
    Phi = y[6]
    # l=2 homogeneous
    Xh = y[7]
    Phi_h = y[8]
    dY = np.zeros_like(y)
    # definte all the intermediate equations that are built into the ODEs

    dY[0] = -(3/5)*(G*MP**(5/3)*M*rho**(1/3))/(K*r**2)
    dY[1] = 4 * np.pi * r**2 * rho
    dY[2] = -rho * G * M / r**2
    # l=0 correction terms
    dY[3] = (2/3)*omg**2*r - (G * M2/r**2)
    dY[4] = 4 * np.pi * r**2 * rho*(dY[0]/dY[2])*p2
    # l=2 particular solution
    dY[5] = -(2*G*M / r**2)*Phi + (8*np.pi/3)*omg**2 * r**3 * G * rho
    dY[6] = ((4*np.pi*r**2*rho/M) - (2/r))*Phi - \
        ((2*X)/(G*M)) + (4*np.pi/(3*M))*rho*omg**2*r**4
    # l=2 homogeneous solution
    dY[7] = -(2*G*M/r**2)*Phi_h
    dY[8] = ((4*np.pi*r**2*rho/M) - (2/r))*Phi_h - ((2*Xh)/(G*M))
    return (dY)


def density_finder(M_f):

    final_mass = []
    rho0 = np.linspace(1e03, 1e06, 200)
    for i in rho0:
        # initial state
        r0 = 0.1  # cm
        M0 = (4/3) * np.pi * r0**3 * i
        P0 = (K/(MP**(5/3))) * i**(5/3)
        # l=0 boundary conditions
        M2_0 = 0
        p2_0 = (1/3) * omg**2 * r0**2
        # l=2 particular boundary conditions
        A = 1
        B = (2*np.pi/3)*G*i*(omg**2 - A)
        Bh = (2*np.pi/3)*G*i*(-A)  # homogeneous bc
        X0 = B * r0**4
        X0_h = Bh * r0**4
        Phi0 = A*r0**2  # same as Phi0_h
        Y0 = np.array([i, M0, P0, p2_0, M2_0, X0, Phi0, X0_h, Phi0])
        t = np.linspace(r0, 5.5e9)

        # integrater
        rk45 = integrate.solve_ivp(
            hydro_ODE, [r0, 5.5e9], Y0, 'RK45', t_eval=t, first_step=0.1)

        mass = rk45.y[1]
        solar_mass = 5e-34 * mass
        corr_mass = 5e-34 * rk45.y[4]
        # get the total mass so central density can be updated
        final_corr_mass = solar_mass + corr_mass
        final_mass.append(final_corr_mass[-1])

    # get the value of rho0 that we need for 0.6M_sun
    real_rho = np.interp(M_f, final_mass, rho0)

    return (final_mass, rho0, real_rho)


########## INSERT FINAL MASS HERE ##############
z = density_finder(0.6)
################################################

print("The Central Density is:", round(z[2], 4))  # central denisty value

# plotting central density vs final mass
plt.figure()
plt.plot(z[1], z[0])
plt.title("Mass v Central Density for Polytropic EOS")
plt.xlabel("Central Density $\\rho_{0}$")
plt.ylabel("Final Mass ($M_{\\odot}$)")
plt.show()


def rk4(rho0):

    # initial state
    r0 = 0.1  # cm
    M0 = (4/3) * np.pi * r0**3 * rho0
    P0 = (K/(MP**(5/3))) * rho0**(5/3)
    M2_0 = 0
    p2_0 = (1/3) * omg**2 * r0**2
    # l=2 particular boundary conditions
    A = 1
    B = (2*np.pi/3)*G*rho0*(omg**2 - A)
    Bh = (2*np.pi/3)*G*rho0*(-A)  # homogeneous bc
    Phi0 = A*r0**2
    X0 = B * r0**4
    X0_h = Bh * r0**4
    Y0 = np.array([rho0, M0, P0, p2_0, M2_0, X0, Phi0, X0_h, Phi0])
    a = 5.5e9  # surface radius
    # a = 1e3
    t = np.linspace(r0, a)  # range of radius values

    # integrator
    rk45 = integrate.solve_ivp(
        hydro_ODE, [r0, a], Y0, 'RK45', t_eval=t, first_step=0.1)

    density = rk45.y[0]
    mass = rk45.y[1]
    solar_mass = 5e-34 * mass
    # print("The Total Solar Mass (uncorrected) is:", solar_mass[-1])
    pressure = rk45.y[2]
    corr_mass = rk45.y[4]
    corr_solar_mass = 5e-34 * corr_mass
    # print("The Corrected Solar Mass term is:", corr_solar_mass[-1])
    corr_pressure = rk45.y[3]
    total_mass = solar_mass + corr_solar_mass
    # print("The Total Mass (corr + uncorr) is:", total_mass[-1])
    total_pressure = pressure + corr_pressure

    # tackle the l=2 solutions - need to match internal and external solutions to solve for constants
    Phi_p = rk45.y[6]
    Phi_h = rk45.y[8]
    X_p = rk45.y[5]
    X_h = rk45.y[7]

    # note: mass is in g not M_sun here
    Aa = np.array([[Phi_h[-1], -1/a**3], [X_h[-1], -G*mass[-1]/(2*a**4)]])
    Bb = np.array([-Phi_p[-1], -X_p[-1]])
    K12 = np.linalg.solve(Aa, Bb)

    # Constants K1 and K2
    K1 = K12[1]
    K2 = K12[0]
    # print("K2=",K2)

    # now assemble phi_in and X_in
    Phi_in = Phi_p + K2 * Phi_h
    X_in = X_p + K2 * X_h

    # calculate zeta 0 and 2
    # array of radii from center to surface of star
    r = np.linspace(r0, a, len(mass))
    zeta0 = (r**2 / (G * mass)) * corr_pressure
    zeta2 = -(r**2 / (G*mass)) * ((1/3)*omg**2*r**2 + Phi_in)

    return (solar_mass, density, r0, pressure, total_mass, corr_pressure, total_pressure, zeta0, zeta2)


qq = rk4(z[2])
# print("Fractional Pressure:",x[3][-1]/x[3][0])

# Now create function to integrate moments of inertia
# one function defining differential equations, one function to integrate
# for 2nd order coharrection, need zeta0 and zeta2, which s mass and p* included in it. Also need d(rho)/dR


def inertia_ODE(r: range, y):
    rho = y[0]
    I0 = y[1]
    M = y[2]
    I2 = y[3]
    dX = np.zeros_like(y)

    # density
    dX[0] = -(3/5)*(G*MP**(5/3)*M*rho**(1/3))/(K*r**2)
    # 0th order moment of inertia
    dX[1] = (8*np.pi/3) * rho * r**4
    # mass
    dX[2] = 4 * np.pi * r**2 * rho
    zeta0, zeta2 = rk4(z[2])[-2:]

    for x in r:
        r.index
    # for i,j in zip(zeta0, zeta2):
    # 2nd order correction moment of inertia
    for i, j in zip(zeta0, zeta2):
        dX[3] = (8*np.pi / 3) * ((1/5)*zeta2[j] - zeta0[i]) * dX[0] * r**4
        return (dX[3])
    return (dX)


def inertia_integrator(rho0):

    # initial state
    r0 = 0.1  # cm
    I0_0 = (8*np.pi/15)*rho0*r0**5
    M0 = (4/3) * np.pi * r0**3 * rho0
    I2_0 = 0
    X0 = np.array([rho0, I0_0, M0, I2_0])
    a = 5.5e9  # surface radius
    # a = 1e3
    t = np.linspace(r0, a)  # range of radius values

    rk45 = integrate.solve_ivp(
        inertia_ODE, [r0, a], X0, 'RK45', t_eval=t, first_step=0.1)

    density = rk45.y[0]
    inertia_0 = rk45.y[1]
    mass = 5e-34 * rk45.y[2]
    inertia_2 = rk45.y[3]
    print("Inertia Integrator final mass:", mass[-1])
    print("0th order MOI = ", inertia_0[-1])
    print("2nd order MOI = ", inertia_2[-1])
    print("total MOI = ", inertia_0[-1]+inertia_2[-1])

    print("Density on surface:", density[-1])

    return (density, inertia_0, mass, inertia_2)


ii = inertia_integrator(z[2])
print(ii[3])


#########################################################################################
################# USE RK4 FUNCTION OUTPUTS FOR ALL MEANINGFUL GRAPHS ####################
#########################################################################################

# Plotting MvR and RhovR for 0th order quantities
t = np.linspace(qq[2], 5.5e9, len(qq[0]))
# t = np.linspace(qq[2], 1e3, len(qq[0]))
''' fig, (ax1, ax2, ax3) = plt.subplots(1,3)
ax1.plot(t, x[0])SeriesData[r, 0, {6.474459949864035*^20 (r^(-1))^Rational[1, 2]
     r^Rational[1, 2]}, -3, 13, 2] ($M^{(0)}_{\\odot})$")
ax2.plot(t, x[1])
ax2.set_title("$\\rho$ vs r")
ax2.set_ylabel('Density (g/cm^3)')
ax2.set_xlabel('Radius (cm)')
ax3.plot(t, x[3])
ax3.set_title("p vs r")
ax3.set_ylabel('Uncorrected Pressure')
ax3.set_xlabel('Radius (cm)')
plt.show()
'''
# plotting total quantities vs R - including second order corrections
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
ax1.plot(t, qq[4])
ax1.set_title("M vs r")
ax1.set_xlabel("Radius (cm)")
ax1.set_ylabel("Total Mass ($M_{\\odot})$")
ax2.plot(t, qq[1])
ax2.set_title("$\\rho$ vs r")
ax2.set_ylabel('Density (g/cm^3)')
ax2.set_xlabel('Radius (cm)')
ax3.plot(t, qq[6])
ax3.set_title("p vs r")
ax3.set_ylabel('Pressure')
ax3.set_xlabel('Radius (cm)')
plt.show()
'''
#0th order moment of inertia v 0 & 2nd order mass
fig, (ax1, ax2) = plt.subplots(1,2)
ax1.plot(qq[0], qq[7]/10**50)
ax1.set_title("$M^{0}$ vs $I^{0}$")
ax1.set_ylabel("I (g*cm^2)")
ax1.set_xlabel("0th Order Mass ($M_{\\odot})$")
ax2.plot(qq[4], qq[7]/10**50)
ax2.set_title("$Total M$ vs $I^{0}$")
ax2.set_xlabel('Total Mass')
ax2.set_ylabel("I (g*cm^2)")
plt.show()
'''
plt.figure()
plt.plot(t, ii[3])
plt.xlabel("Radius (cm)")
plt.ylabel("2nd order MoI (g*cm^2)")
plt.title("I2 vs Radius of Star for Small R")
plt.show()
