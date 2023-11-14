import numpy as np
import math
from priors import prior_draw
from HelperCalculations import derivative
from graphing import makeCornerPlot

#test case - resample priors for A and Phi
A = 0.5
phi = np.pi / 3.
params = [A, phi]
x_real = A*np.cos(phi)
y_real = A*np.sin(phi)
corner_params = [x_real, y_real]
#prior limits
A_low_lims = 0
A_high_lims = np.sqrt(2)

phi_low_lims = -np.pi
phi_high_lims = np.pi

low_lims = [A_low_lims, phi_low_lims]
high_lims = [A_high_lims, phi_high_lims]

def testResampling(params, nsteps):
    A_prior = []
    phi_prior = []
    i=0
    while i < nsteps:
        prior = prior_draw(high_lims, low_lims, params)
        A_prior.append(prior[0])
        phi_prior.append(prior[1])
        i += 1
    
    #two cases - 1. jacobian weight in corner plot. 2. eqn (17): p(x,y) = (1/A)p(A)p(phi)
    #now convert to x and y
    x_prior = []
    y_prior = []
    jacobian = []
    AA = []
    for i in range(len(A_prior)):
        x = A_prior[i] * np.cos(phi_prior[i])
        y = A_prior[i] * np.sin(phi_prior[i])
        if np.abs(x) < 1 and np.abs(y) < 1:
            x_prior.append(x)
            y_prior.append(y)
            j = get_Jacobian_test([x, y])
            jacobian.append(j)
            AA.append((x**2 + y**2)**(0.5))
    makeCornerPlot(A_prior, phi_prior, params, "A", r"$\phi$")
    makeCornerPlot(x_prior, y_prior, corner_params, "x", "y")
    makeCornerPlot(x_prior, y_prior, corner_params, "x", "y", weight=AA)

    return

def get_Jacobian_test(p):
    step_size = [1.e-10, 1.e-10]
    jacobian = np.zeros((np.size(p), np.size(p)))
    for i in range(len(p)):
        for j in range(len(p)):
            fx = lambda params: testParameterization(params)[i]
            dF = lambda x, params: derivative(fx, params, step_size, x)
            jacobian[i, j] = dF(j, p)
    #invert this matrix, need d(m1m2)/d(b,d)
    inverseJac = np.linalg.inv(jacobian)
    #take determinant of jacobian
    detjac = np.abs(np.linalg.det(inverseJac))
    #detjac = np.abs(np.linalg.det(jacobian))
    return detjac

def testParameterization(params):
    x = params[0]
    y = params[1]
    A = (x**2 + y**2)**(0.5)
    phi = np.arctan(y / x)
    return(A, phi)

testResampling(params, 100000)