import numpy as np
from Calculations import Frequency_1PN
from const import *

def main(freq0, chirpMass, totalMass):
    #Calculates frequency derivatives and other relevant quantities for binary white dwarf in LISA

    fdot_1PN, fddot_1PN = Frequency_1PN(freq0, chirpMass, totalMass)
    print("1PN Freq Derivative =", fdot_1PN)
    print("1PN Freq Second Derivative =", fddot_1PN)

main(10.e-3, 0.522*MSOLAR, 1.4*MSOLAR)