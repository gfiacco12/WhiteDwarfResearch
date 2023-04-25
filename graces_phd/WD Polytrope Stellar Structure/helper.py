def getScaledTotalMass(M, M2, omega_new, omega_old):
        return ((((omega_new/omega_old)**2) * M2) + M) * 5e-34

def getScaledTotalInertia(I0, I2, omega_new, omega_old):
        return ((((omega_new/omega_old)**2) * I2) + I0)