#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Rafael Figueroa
Hybrid Automata Model and simulation using HaSimPy
"""

from robot_vars import *
import numpy as np
from hasimpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from feature_calculations import psiMatrix

sin = np.sin
cos = np.cos
tan = np.tan
pi = np.pi


def w0(X):
    featMatrix = f0(X)
    matrixTrans = np.transpose(featMatrix)
    wVal = (1/(matrixTrans * featMatrix) + (10-6)) * matrixTrans
    return wVal


def hw(X,z):
    zw = (r0(X)) * w0(X)
    zPsi =  zw * f0(X)
    hwSum = sum(zPsi)
    return hwSum


def tractor_system(g,y,w,z):
    linCont = ((g-y) - y + hw()) + (sum(((z)*w))/sum((z)))
    return linCont


def euler_integration(f, a, a2, b, n):
    c = (a2 - a)/float(n)
    x = a
    y = b
    for i in range(n):
        y += c * f(x,y)
        x += c
    return y
# http://code.activestate.com/recipes/577647-ode-solver-using-euler-method/

lengthArray = [3, 5, 6]
loopNum = [0, 1, 2]
limbAngleArray = [2, 4, 6]
prevAngleCos = 0
prevAngleSin = 0


def p1(lengthArray, loopNum, limbAngleArray, prevAngleCos, prevAngleSin):
    Length = ()
    for i in loopNum:
        angleCos = cos(limbAngleArray[i])
        angleSin = sin(limbAngleArray[i])
        angleArray = np.array([(angleCos + prevAngleCos), (angleSin + prevAngleSin)])
        PArray = Length * angleArray
        prevAngleCos = angleCos
        prevAngleSin = angleSin
        return PArray

timeStamps = [0, 1, 2, 3, 4]


def psiMatrix(k, T, psiMatrix):
    kcount = len(k)
    timecount = len(T)
    psiVal = [[0 for x in range(kcount)] for y in range(timecount)]
    for i in kcount:
        for j in timecount:
            psiVal[i][j] = np.array()

    return psiVal


