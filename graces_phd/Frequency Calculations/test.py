import numpy as np
import math

chirpMass = [1, 1, 1, float('nan'), 1]
totalMass = [2, -1, -0.5, 2, 2]
print(chirpMass)
print(totalMass)
chirpMass_new = []
totalMass_new = []

for mass in range(len(chirpMass)):
    if not math.isnan(chirpMass[mass]) and totalMass[mass] > 0:
        chirpMass_new.append(chirpMass[mass])
        totalMass_new.append(totalMass[mass])
print(chirpMass_new)
print(totalMass_new)