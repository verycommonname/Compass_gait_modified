from robot_vars import *
import numpy as np
from hasimpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from phase_calculation import zValue

c = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]

exp = np.exp

def featureCalc(zValue, c, h):
    countNum = [None] * len(zValue)
    calcZ = [countNum]
    psiNum = [None] * len(psiVal)
    psiTotal = [None] * len(psiVal)
    for i in countNum:
        calcZ.append(exp(-0.5((zValue[i] - c[i])**2)/ h))
        #calcZ[i] = np.array(exp(-0.5((z - c[i])**2)/ h))

        plt.plot(psiTotal)
        plt.ylabel('psi value')
        plt.xlabel('z value')

    for j in psiNum:
        psiTotal.append(calcZ[j]/sum(calcZ[j]))*j
        #psiTotal[j] = np.array(calcZ[j]/sum(calcZ[j]))*j

    plt.show()

    return psiTotal

psiNorm = featureCalc(zValue, c, h)


def psiMatrix(k, T, psiNorm):
    kcount = [None] * len(k)
    timecount = [None] * len(T)
    psiVal = [[0 for x in range(kcount)] for y in range(timecount)]
    for i in kcount:
        for j in timecount:
            psiVal[i][j].append(psiNorm[i])
            #psiVal[i][j] = np.array(psiNorm[i])

    return psiValue

psiMtx = psiMatrix(k, T, psiNorm)

