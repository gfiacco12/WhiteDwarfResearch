from numpy import number
from calculations import findDensity, integrateStar, integrateInertia


def main(r0: number, a: number):
    # Function to calculate central density for Polytropic White Dwarf to achieve specific final mass
    # INPUTS: Desired final mass (in Solar Mass Units), range of central density values (start, stop, number of steps)
    # INPUTS CONT: Initial radius r0, surface radius of star
    rho = findDensity(0.6, 1e03, 1e06, 200, r0, a)
    print("The Central Density is:", round(rho, 4))  # central denisty value

    # integrates the stellar structure for star using proper central density
    K2 = integrateStar(rho, r0, a)

    integrateInertia(rho, r0, a, K2)

main(0.1, 5.5e9)
