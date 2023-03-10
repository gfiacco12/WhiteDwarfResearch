from abc import abstractmethod
from numpy import number
import numpy as np
from scipy import integrate

# class definition: "custom type" -- structure made up of different types(int, string, num, etc). Lets you have
# multiple different values on a single object. Lets you run functions as well.
# Class is like a template, each time you call a class(aka create an instance), it creates a copy...
# that has its own specific variables

# all class functions start with self.


class State:
    # constructor -- can call when creating an instance. Used to set default values if needed
    # called whenever you call class as a function, used to create instance
    def __init__(self, rho: number = 0):
        self.rho: number = rho
        self.stepIndex: number = 0

    def integrateSelf(self, cls: type, start: number, end: number, stop_condition=None):
        self.stepIndex = 0
        
        t = np.linspace(start, end)

        paramNames = self.getIntegrationParamNames()
        paramValues = self.getIntegrationParamValues()
        y0 = np.array(paramValues)

        #impliments a version of the integrator using stop conditions, so it stops when p < 0
        #instead of when we hit a specified radius
        rk45 = None
        if (stop_condition == None):
            rk45 = integrate.solve_ivp(
                self.ODE, [start, end], y0, t_eval=t, first_step=start)
        else:
            rk45 = integrate.solve_ivp(
                self.ODE, [start, end], y0, t_eval=t, first_step=start, atol=1e-8, events=stop_condition)

        # creates a list of states
        length = len(rk45.y[0])
        states: list[cls] = []  # creates empty list
        # for every step the integrator takes, we create an object with all the values
        # makes an array by step, not by parameter, and puts each array in list state
        for i in range(length):
            stateEntry = cls()
            j = 0
            for param in paramNames:
                setattr(stateEntry, param, rk45.y[j][i])
                j = j+1
            states.append(stateEntry)

        #null checks to make sure everyhing exists to get the radius on termination
        if (rk45.t_events != None and len(rk45.t_events) > 0 and len(rk45.t_events[0]) > 0):
            return states, rk45.t_events[0][0]
        
        return states

    @abstractmethod
    # initial state for all the differential equations (starting values at r0)
    def setInitialState(self, r0: number):
        raise NotImplementedError()

    @abstractmethod
    def getIntegrationParamNames(self):
        raise NotImplementedError()

    @abstractmethod
    def getIntegrationParamValues(self):
        raise NotImplementedError()

    @abstractmethod
    def ODE(self, r: number, y):
        raise NotImplementedError()