import numpy as np
from const import MSOLAR, SECSYEAR
from priors import prior_draw
from HelperCalculations import derivative, getChirpMass
from graphing import makeCornerPlot
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.signal import windows

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

#test FFT stuff:
def FFT_test():
    #h(t)1
    f0 = 20.e-3
    t_obs = 0.5 * SECSYEAR
    t = np.arange(0, t_obs, 2**4)
    h_t = np.sin(2*np.pi*f0*t)

    #h2(t)
    month = (1/12)*SECSYEAR
    df = 1 / (month**2)
    print(df)
    h2t = np.sin(2*np.pi*((f0*t)+(df*(t**2))))

    N = len(h_t)
    n_half = int(N/2)
    h_f = fft(h_t)
    f = fftfreq(N, 2**4)
    win = windows.tukey(N, 0.1)
    h_half = h_f[:n_half]
    f_half = f[:n_half]

    h2f = fft(h2t) * win
    h2f_half = h2f[:n_half]

    plt.figure()
    plt.loglog(f_half, f_half*np.abs(h_half), label="h(f)")
    plt.loglog(f_half, f_half*np.abs(h2f_half), label="h2(f)")
    plt.legend()
    plt.xlabel("Hz")
    plt.ylabel("h(f)")
    plt.xlim(19.99e-3, 20.01e-3)
    #plt.axvline(x = 0.02, color = 'r')
    plt.show()
    return
FFT_test()