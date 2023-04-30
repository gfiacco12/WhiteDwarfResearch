def getScaledMass2(M2, omega_new, omega_old):
        return (((omega_new/omega_old)**2) * M2)

def getScaledTotalMass(M1, M2):
        return (M1 + M2)

def getScaledInertia2(I2, omega_new, omega_old):
        return (((omega_new/omega_old)**2) * I2) 

def getScaledTotalInertia(I0, I2):
        return (I0 + I2)