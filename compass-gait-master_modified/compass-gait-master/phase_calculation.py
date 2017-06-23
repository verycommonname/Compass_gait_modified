from robot_vars import *
import numpy as np
from hasimpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation

timeStamps = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

alpha = 0.1

#plt.ylim(0, 1)

def phaseCalc(timeStamps, alpha):
    zVal = [None] * len(timeStamps)
    for i in range(len(timeStamps)):
        j = timeStamps[i] - 1
        if j is 0:
            zVal.append(1) #zVal[j] = 1
        else:
            alphaVal = np.gradient(-alpha) * timeStamps[i]
            timeVal = np.gradient(timeStamps[i]) * -alpha
            zVal.append(timeStamps[i] + (alphaVal + timeVal))
            # zVal.append(timeStamps[i] + np.gradient(-alpha*timeStamps[i]))
            # (timeStamps[i])# + np.diff(-alpha*timeStamps[i]))
            #zVal[j] = np.array(timeStamps[i] + np.diff(-alpha*timeStamps[i]))

    return zVal

zValue = phaseCalc(timeStamps, alpha)

plt.plot(zValue)
plt.ylabel('z value')
plt.xlabel('t value')
plt.show()