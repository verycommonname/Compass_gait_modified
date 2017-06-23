from robot_vars import *
import numpy as np
from hasimpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from feature_calculations import psiMtx


def linAlg(psiMtx, alpha, I, F):
    aVal = ((np.transpose(psiMtx)*psiMtx)+(alpha*I))**-1 # add power of minus one
    bVal = (np.transpose(psiMtx)) * F
    linAlg = np.linalg.solve(aVal,bVal)
    return linAlg

linearAlg = linAlg(psiMtx, alpha, I, F)

plt.plot(linearAlg)
plt.ylabel('I value')
plt.xlabel('t value')
plt.show()