from graphing import plotCentralDensity, plotStellarStructure
from states.state import State
from states.hydroState import HydroState
from states.inertiaState import InertiaState
from numpy import number
import numpy as np
from globals import *

# Density finder function

def findDensity(finalSolarMass: number, start: number, end: number, steps: number, r0: number, a: number):
    final_masses = []  # list to store all final solar masses for each step in rho0
    rho0 = np.linspace(start, end, steps)
    for rho in rho0:
        # makes new state and lets you put in variable "rho"
        state = HydroState(rho)
        # calling function setInitialState in the instance created above
        state.setInitialState(r0)
        # calls integration function
        integrationResults: list[HydroState] = state.integrateSelf(HydroState, r0, a)
        # adds specifically last value for each final_corr_mass array
        final_masses.append(integrationResults[-1].TotalMass)

    # get the value of rho0 that we need for 0.6M_sun
    real_rho = np.interp(finalSolarMass, final_masses, rho0)

    plotCentralDensity(final_masses, rho0)

    return real_rho


def integrateStar(rho: number, r0: number, a: number):
    state = HydroState(rho)
    state.setInitialState(r0)
    integrationResults: list[HydroState] = state.integrateSelf(HydroState, r0, a)
    finalStep = integrationResults[-1]

    # note: mass is in g not M_sun here
    Aa = np.array([[finalStep.Phi_h, -1/a**3],
                  [finalStep.X_h, -G*finalStep.M/(2*a**4)]])
    Bb = np.array([-finalStep.Phi_p, -finalStep.X_p])
    K12 = np.linalg.solve(Aa, Bb)

    # Constants K1 and K2
    K1 = K12[1]
    K2 = K12[0]

    # calculate zeta 0 and 2
    # array of radii from center to surface of star
    rArr = np.linspace(r0, a, len(integrationResults))
    zeta0 = []
    zeta2 = []
    densities = []
    totalMasses = []
    totalPressures = []
    for i in range(len(integrationResults)):
        step = integrationResults[i]
        Phi_in = step.Phi_p + K2 * step.Phi_h
        # X_in = step.X_p + K2 * step.X_h
        r = rArr[i]
        z0 = (r**2 / (G * step.M)) * step.P2
        z2 = -(r**2 / (G*step.M)) * ((1/3)*omg**2*r**2 + Phi_in)
        zeta0.append(z0)
        zeta2.append(z2)
        densities.append(step.rho)
        totalMasses.append(step.TotalMass)
        totalPressures.append(step.TotalPressure)

    t = np.linspace(r0, a, len(integrationResults))
    
    plotStellarStructure(t, densities, totalMasses, totalPressures)

    return K2

def integrateInertia(rho: number, r0: number, a: number, k2: number):
    state = InertiaState(rho)
    state.setInitialState(r0, k2)
    integrationResults: list[InertiaState] = state.integrateSelf(InertiaState, r0, a)
    finalStep = integrationResults[-1]

    print("Inertia Integrator final mass:", finalStep.SolarMass)
    print("0th order MOI = ", finalStep.I0)
    print("2nd order MOI = ", finalStep.I2)
    print("total MOI = ", finalStep.TotalI)