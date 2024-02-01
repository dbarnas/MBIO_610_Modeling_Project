## IMPORT LIBRARIES
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import rc
import seaborn as sns
import pandas as pd
from scipy import stats
from scipy.optimize import curve_fit
from scipy.integrate import odeint
from scipy.integrate import solve_ivp

# My model:
# Set parameter values for testing (using constant sgd impact factor)
rC = 0.2  # coral growth
rM = 0.4  # algae growth

aC = 0.2  # algal impact onto coral
aM = 0.2  # coral impact onto algae

nC = 1.5  # coral impact from sgd (dummy factors for phosphate impact)
nM = 0.5  # algae impact from sgd

K = 1.0  # max carrying capacity for coral and algae in the absence of second group

dt = 0.1  # step size
NUMYEARS = 100
NUMSTEPS = int(NUMYEARS / dt)

# Create a dictionary object called 'parameters_dict'
parameters_dict = {'dt': dt,
                   'NUMSTEPS': NUMSTEPS,
                   'rC': rC,
                   'rM': rM,
                   'aC': aC,
                   'aM': aM,
                   'nC': nC,
                   'nM': nM,
                   'K': K
                   }

# set initial conditions for C0 and M0
C0 = 0.5
M0 = 0.5

# My model:
def dNdt(C, M, P):
    dt = P['dt']

    # Calculate the derivative
    dC = (rC * C / K) * (K - C - aC * M - K * nC * C) * dt  # Coral equation
    dM = (rM * M / K) * (K - M - aM * C - K * nM * M) * dt  # Macroalgae equation

    return dC, dM


def RK2(C, M, P):  # 2nd-order Runge-Kutta

    C_init = C
    M_init = M

    dC1, dM1 = dNdt(C, M, P)

    C1 = C + 0.5 * dC1  # what is the 0.5 value from??
    M1 = M + 0.5 * dM1

    dC2, dM2 = dNdt(C1, M1, P)

    dCave = (dC1 + dC2) / 2
    dMave = (dM1 + dM2) / 2

    C = C_init + dCave
    M = M_init + dMave

    return C, M


def run_model_RK2(INIT_C, INIT_M, P):
    NUMSTEPS = P['NUMSTEPS']

    C = np.zeros((NUMSTEPS))
    M = np.zeros((NUMSTEPS))

    C[0] = INIT_C
    M[0] = INIT_M

    for step in np.arange(0, NUMSTEPS - 1):
        C[step + 1], M[step + 1] = RK2(C[step], M[step], P)

    return C, M



def main():
    global nC, nM
    # Generate output

    # for loop to create random initial population sizes for both C and M (0 < N < 1)
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 6))

    c1 = 'yellow'
    c2 = 'red'
    n = 102

    cmap = plt.get_cmap()

    for i in np.arange(0, 100):
        # randomize the sgd parameter, while holding starting population size constant
        nC = np.random.random()  # populates values b/n 0 and 1 - may need to revisit to get larger or neg values
        nM = np.random.random()

        C_array, M_array = run_model_RK2(C0, M0, parameters_dict)

        ax1.plot(C_array, c=cmap(nC), alpha=0.5)  # plot the coral time series in plot C
        #     ax2.plot(M_array, c='green', alpha = 0.5) # plot the macroalgae time series in plot C
        #    ax1.plot(C_array, color=get_color_gradient("#fffc00", "#ff0300", 100), alpha = 0.5) # plot the coral time series in plot C
        ax2.plot(M_array, c=cmap(nM), alpha=0.5)  # plot the macroalgae time series in plot D

    # ax1.legend(['SGD impact \nfactor range'], loc='lower left', bbox_to_anchor=(-0.75, 0.8), fontsize=14)
    ax1.tick_params(axis='x', labelsize=14)
    ax1.tick_params(axis='y', labelsize=14)
    ax1.set_xlabel("Time", fontsize=14)
    ax1.set_ylabel("Relative abundance", fontsize=14)
    ax1.set_title('[Coral] SGD impact on \ncoral reef community dynamics', fontsize=14)
    # Add colorbar, make sure to specify tick locations to match desired ticklabels

    fig.colorbar(mpl.cm.ScalarMappable(norm=None, cmap=cmap), ax=ax2)
    #ax2.legend(['SGD impact \nfactor range'], loc='lower left', fontsize=14)
    ax2.tick_params(axis='x', labelsize=14)
    ax2.tick_params(axis='y', labelsize=14)
    ax2.set_xlabel("Time", fontsize=14)
    # ax2.set_ylabel("Relative abundance", fontsize=14);
    ax2.set_title('[Macroalgae+Turf] SGD impact on \ncoral reef community dynamics', fontsize=14)
    plt.show()


if __name__ == "__main__":
    main()
