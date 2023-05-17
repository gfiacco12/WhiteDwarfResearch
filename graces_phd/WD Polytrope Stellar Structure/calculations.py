from helper import getScaledMass2, getScaledInertia2, getScaledTotalMass, getScaledTotalInertia
from graphing import plotCentralDensity, plotMassRadiusRelation, plotStellarStructure, plotMassInertiaRelation, plot3DMassRadiusVelocity,plot3DMassInertiaVelocity,plot2DMassInertiaVelocity
from states.state import State
from states.hydroState import HydroState
from states.inertiaState import InertiaState
from numpy import number
import numpy as np
from globals import *
import matplotlib.pyplot as plt

# Density finder function

#implimenting the stop conditions -- looks for sign changes, not 0 to terminate
def eventAttr():
    def decorator(func):
        func.direction = 0
        func.terminal = True
        return func
    return decorator

@eventAttr()
def stop_condition(r, y):
    if (y[0] < 100):
        return -1
    return y[0]

def findDensity(finalSolarMass: number, start: number, end: number, steps: number, r0: number, a: number):
    final_r_withStop = []
    final_masses_withStop = []
    final_masses = []  # list to store all final solar masses for each step in rho0
    final_I_withStop = []

    rho0 = np.linspace(start, end, steps)
    for rhoI in rho0:
        # makes new state and lets you put in variable "rho"
        state = HydroState(rhoI)
        # calling function setInitialState in the instance created above
        state.setInitialState(r0)
        # calls integration function
        integrationResults: list[HydroState] = state.integrateSelf(HydroState, r0, a)

        # adds specifically last value for each final_corr_mass array
        final_masses.append(integrationResults[-1].TotalMass)

        #trying to un-fix the radius, so mass and radius can vary with central
        #density to get MvsR
        #also get 0th order moment of inertia
        integrationResultsWithStop: list[HydroState]
        t: number
        integrationResultsWithStop, t = state.integrateSelf(HydroState, r0, a, stop_condition)
        final_masses_withStop.append(integrationResultsWithStop[-1].TotalMass)
        final_r_withStop.append(t)
        final_I_withStop.append(integrationResultsWithStop[-1].I0)

    # get the value of rho0 that we need for 0.6M_sun
    real_rho = np.interp(finalSolarMass, final_masses, rho0)

   # plotCentralDensity(final_masses, rho0)
    #plotMassRadiusRelation(final_masses_withStop, final_r_withStop)
    # plotMassInertiaRelation(final_masses_withStop, final_I_withStop)

    return real_rho

def integrateStar(rho: number, r0: number, a: number):
    initialstate = HydroState(rho)
    initialstate.setInitialState(r0)
    integrationResults: list[HydroState] = initialstate.integrateSelf(HydroState, r0, a)
    t = np.linspace(r0, a, len(integrationResults))
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

    print("M0 = ",finalStep.M* 5e-34)
    print("M2 = ", finalStep.M2* 5e-34)
    
    #plotStellarStructure(t, densities, totalMasses, totalPressures)

    return K2

def integrateInertia(rho: number, r0: number, a: number, k2: number):
    initialstate = InertiaState(rho)
    initialstate.setInitialState(r0, k2)
    integrationResults: list[InertiaState] = initialstate.integrateSelf(InertiaState, r0, a)
    finalStep = integrationResults[-1]

    #print("Inertia Integrator final mass:", finalStep.SolarMass)
    print("0th order MOI = ", finalStep.I0)
    print("2nd order MOI = ", finalStep.I2)
    print("total MOI = ", finalStep.TotalI)

    t = np.linspace(r0, a, len(integrationResults))
    I_0 = []
    I_2 = []
    I_total = []
    for i in range(len(integrationResults)):
        step = integrationResults[i]
        I_0.append(step.I0)
        I_2.append(step.I2)
        I_total.append(step.TotalI)


def getFinalMassesVsFinalInertia(start: number, end: number, steps: number, r0: number, a: number, K2: number):
    final_masses_withStop = []
    final_m_withStop = []
    final_m2_withStop = []
    final_I_withStop = []
    final_I0_withStop = []
    final_I2_withStop = []

    rho0 = np.linspace(start, end, steps)
    for rhoI in rho0:
        # makes new state and lets you put in variable "rho"
        initialstate = InertiaState(rhoI)
        # calling function setInitialState in the instance created above
        initialstate.setInitialState(r0, K2)

        #trying to un-fix the radius, so mass and radius can vary with central
        #density to get MvsR
        #also get 0th order moment of inertia
        integrationResultsWithStop: list[InertiaState]
        integrationResultsWithStop, t = initialstate.integrateSelf(InertiaState, r0, a, stop_condition)
        final_masses_withStop.append(integrationResultsWithStop[-1].TotalMass)
        final_m_withStop.append(integrationResultsWithStop[-1].SolarMass)
        final_m2_withStop.append(integrationResultsWithStop[-1].SolarMass2)
        final_I_withStop.append(integrationResultsWithStop[-1].TotalI)
        final_I0_withStop.append(integrationResultsWithStop[-1].I0)
        final_I2_withStop.append(integrationResultsWithStop[-1].I2)
    
    #obtain scaling for angular velocity for M and I2
    #omega = np.linspace(omg, 0, 3) # Note: eventually convert this to frequencies not angular
    omega = np.linspace(omg, 0, 11)
    finalScaledMasses = []
    finalScaledM2 = []
    finalScaledInertia = []
    finalScaledI2 = []

    for i in range(len(final_masses_withStop)):
        finalScaledMasses.append([])
        finalScaledInertia.append([])
        finalScaledM2.append([])
        finalScaledI2.append([])
        M2_scaled = final_m2_withStop[i]
        I2_scaled = final_I2_withStop[i]
        for j in range(len(omega)):
            if j != 0:
                omega_j = omega[j]
                omega_prev = omega[j-1]
                M2_scaled = getScaledMass2(M2_scaled, omega_j, omega_prev)
                I2_scaled = getScaledInertia2(I2_scaled, omega_j, omega_prev)
            finalScaledM2[i].append(M2_scaled)
            scaledMasses = getScaledTotalMass(final_m_withStop[i], M2_scaled)
            finalScaledMasses[i].append(scaledMasses)
            finalScaledI2[i].append(I2_scaled)
            scaledInertia = getScaledTotalInertia(final_I0_withStop[i], I2_scaled)
            finalScaledInertia[i].append(scaledInertia)

    plot2DMassInertiaVelocity(finalScaledMasses, finalScaledI2)

    #plotMassInertiaRelation(final_masses_withStop, final_I_withStop)
