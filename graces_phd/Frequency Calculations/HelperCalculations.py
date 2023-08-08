from numpy import number
import numpy as np
from const import *
import scipy as sc
from scipy import integrate, linalg
import matplotlib.pyplot as plt

def getChirpMass(mass1, mass2):
    return((mass1*mass2)**(3/5) / (mass1 + mass2)**(1/5))

def getTotalMass(mass1, mass2):
    return(mass1 + mass2)

