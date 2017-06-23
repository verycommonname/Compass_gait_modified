from robot_vars import *
import numpy as np
from hasimpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
# import psi val from additional calculations


def targetCalc(z, w, psiVal, alpha, beta, g, y, f, t):
    psiSum = [None] * len(psiVal)
    tSum = [None] * len(t)
    zFunct = [None]
    tFunct = [None]
    for i in psiSum:
        zFunct.append(sum(psiVal[i] * w[i]))

    for j in tSum:
        velosity[j].append((y * (j + 1)) - (y*j))/dj
        trajectory[j] = (alpha * (beta + (g - y) - velosity[j])) + f
        tFunct[j].append((trajectory[j] - (alpha(beta(g - y[j])))) - velocity[j])

    return
