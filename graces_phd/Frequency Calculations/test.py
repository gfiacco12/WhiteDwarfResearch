import numpy as np
import math
from const import MSOLAR
from priors import prior_draw
from HelperCalculations import derivative, getChirpMass
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

#testResampling(params, 100000)

def test_mass_resampling(params, nsteps):
    m1_prior = []
    q_prior = []
    high_lims_mass = [1.4, 1.0]
    low_lims_mass = [0.3, 0.]
    i=0
    while i < nsteps:
        prior = prior_draw(high_lims_mass, low_lims_mass, params)
        m1_prior.append(prior[0])
        q_prior.append(prior[1])
        i += 1
    #apply mass cuts
    m1 = []
    q = []
    true_chirp = getChirpMass(0.7, 0.6)
    print(true_chirp)
    print((0.857143**(3/5)*0.7) /(1 + 0.857143)**(1/5) )
    for j in range(len(m1_prior)):      
        m2 = q_prior[j] * m1_prior[j]
        mc = (m1_prior[j] * q_prior[j]**(3./5.)) / (1 + q_prior[j])**(1./5.)
        if m2 >= 0.3 and m2 <= m1_prior[j]:
            if mc <= (true_chirp + 0.005*true_chirp) and mc >= (true_chirp - 0.005*true_chirp):
                m1.append(m1_prior[j])
                q.append(q_prior[j])

    masses_Msun = [0.7, 0.857143]
    makeCornerPlot([m1, q], masses_Msun, labels=[r"$M_{1}$", "q"])
    return
mass_param = [0.7, 0.857143]
#test_mass_resampling(mass_param, 500000)

print(getChirpMass(0.3, 0.3))
print(getChirpMass(1.4, 1.4))