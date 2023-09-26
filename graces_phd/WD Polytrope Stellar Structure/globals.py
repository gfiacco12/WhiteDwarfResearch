import numpy as np

# start by defining global variables needed for all the equations - in cgs ick
# Constants needed: G, mp, mu, kb, sig_sb,hbar, electron mass
G = 6.67e-8
MP = 1.6726e-24  # g
MU = 0.8
KB = 1.38e-16
SIG = 5.67e-5
H = 6.626e-27
ME = 9.1094e-28
K = (H ** 2 / (20*ME)) * (3/(np.pi))**(2/3) * (1/2)**(5/3)
#K = (H ** 2 / (5*ME)) * (3/(8*np.pi))**(2/3)

re = 5e9 + 41123358.4 - (1/2)*(-1682956259) #equatorial radius for 0.6M
#omg = (G*1.18937e33/((5e9)**3))**(1/2)
omg = np.pi * 9.e-3
print(omg/np.pi)
omg_K = (G*1.193e33/((re)**3))**(1/2)  # angular velocity - keplerian
