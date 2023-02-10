from typing import Callable
from states.state import State
from numpy import number
from globals import *

# class definition: "custom type" -- structure made up of different types(int, string, num, etc). Lets you have
# multiple different values on a single object. Lets you run functions as well.
# Class is like a template, each time you call a class(aka create an instance), it creates a copy...
# that has its own specific variables

# all class functions start with self.


class HydroState(State):
    @property
    def TotalPressure(self) -> number:
        return self.P + self.P2

    @property
    def TotalMass(self) -> number:
        return (self.M + self.M2) * 5e-34

    # constructor -- can call when creating an instance. Used to set default values if needed
    # called whenever you call class as a function, used to create instance
    def __init__(self, rho: number = 0):
        super().__init__(rho)

    # initial state for all the differential equations (starting values at r0)
    def setInitialState(self, r0: number):
        rho = self.rho
        # initial state
        self.M: number = (4/3) * np.pi * r0**3 * rho
        self.P: number = (K/(MP**(5/3))) * rho**(5/3)
        # l=0 boundary conditions
        self.M2: number = 0
        self.P2: number = (1/3) * omg**2 * r0**2
        # l=2 particular boundary conditions
        A = 1
        B = (2*np.pi/3)*G*rho*(omg**2 - A)
        Bh = (2*np.pi/3)*G*rho*(-A)  # homogeneous bc
        self.X_p: number = B * r0**4
        self.X_h: number = Bh * r0**4
        self.Phi_h: number = A*r0**2
        self.Phi_p: number = self.Phi_h

    def getIntegrationParamNames(self):
        return ["rho", "M", "P", "P2", "M2", "X_p", "Phi_p", "X_h", "Phi_h"]

    def getIntegrationParamValues(self):
        return [self.rho, self.M, self.P, self.P2, self.M2, self.X_p, self.Phi_p, self.X_h, self.Phi_h]

    def ODE(self, r: number, y):
        dY = np.zeros_like(y)

        rho = y[0]
        M = y[1]
        # l=0
        p2 = y[3]
        M2 = y[4]
        # l=2, particular
        X = y[5]
        Phi = y[6]
        # l=2 homogeneous
        Xh = y[7]
        Phi_h = y[8]

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