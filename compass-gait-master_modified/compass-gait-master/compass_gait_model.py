#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Rafael Figueroa
Hybrid Automata Model and simulation using HaSimPy
"""

from robot_vars import *
import numpy as np
from hasimpy import *

sin = np.sin
cos = np.cos
tan = np.tan
pi = np.pi


def u0(X = None):
    """Passive Walking"""
    return 0.0


def f0(X, u = u0(), t=0):
    """Compass Gait Dynamics"""

    hsw = X[0]
    hst = X[1]
    hdsw = X[2]
    hdst = X[3]
    u = u()

    dx1 = hdsw
    dx2 = hdst

    dx3 = (b*l*m*(b*hdsw**2*l*m*sin(hst - hsw) + \
           g*(M*l + a*m + l*m)*sin(hst) - 1.0*u)*cos(hst - hsw) - \
           (a**2*m + l**2*(M + m))*(b*g*m*sin(hsw) + \
           b*hdst**2*l*m*sin(hst - hsw) - 1.0*u))/(b**2*m*(a**2*m - \
           l**2*m*cos(hst - hsw)**2 + l**2*(M + m)))

    dx4 = (b*(b*hdsw**2*l*m*sin(hst - hsw) + g*(M*l + a*m + l*m)*sin(hst) - \
          1.0*u) - l*(b*g*m*sin(hsw) + b*hdst**2*l*m*sin(hst - hsw) - \
          1.0*u)*cos(hst - hsw))/(b*(a**2*m - l**2*m*cos(hst - hsw)**2 +  \
          l**2*(M + m)))

    dX = np.array([dx1, dx2, dx3, dx4])

    return dX 


def r0(X):
    """Reset Map at impact"""

    hsw = X[0]
    hst = X[1]
    hdsw = X[2]
    hdst = X[3]

    print ('reset map Xold:', X)

    Pl = [[-a*l*m*cos(1.0*hst - hsw)/(M*l**2 + a**2*m - l**2*m*cos(1.0*hst - \
        hsw)**2 + l**2*m), (-M*a*l**2 + M*l**3*cos(1.0*hst - hsw)**2 - \
        a**3*m + 2*a*l**2*m*cos(1.0*hst - hsw)**2 - a*l**2*m)/(b*(M*l**2 + \
        a**2*m - l**2*m*cos(1.0*hst - hsw)**2 + l**2*m))], \
        [-a*b*m/(M*l**2 + a**2*m - l**2*m*cos(1.0*hst - hsw)**2 + l**2*m), \
        l*(M*l + a*m)*cos(1.0*hst - hsw)/(M*l**2 + a**2*m - \
        l**2*m*cos(1.0*hst - hsw)**2 + l**2*m)]]

    P = np.array(Pl)
    J = np.array([[0, 1], [1, 0]])

    print ('reset map P:', P)

    # q = [hsw, hst] transposed
    # dq = [hdsw, hdst] transposed
    # q = J*q (before and after)
    # dq = P*dq (before and after)
    
    q_old = np.array([hsw, hst])
    dq_old = np.array([hdsw, hdst])

    q = np.dot(J, q_old.transpose())
    dq = np.dot(P, dq_old.transpose())

    Xnew = np.hstack((q, dq))
    print ('reset map Xnew:', Xnew)

    return Xnew


def g0(X):
    """Guard: activates (True) when the swing leg impacts the ground"""
    hsw = X[0]
    hst = X[1]
   
    impact_cond = tolEqual(tan(hst/2 + hsw/2), tan(pi-gamma))
    walking_cond = (hsw > hst)
    print ('hs',tan(hst/2 + hsw/2), 'gm', tan(pi-gamma), 'imp', impact_cond, 'wal', walking_cond)
    
    return impact_cond and walking_cond


def avoid(X):
    """No avoid set analysis for this system"""
    return False  # Never goes into the avoid set





# ---- Animation and plots ---- #
e0 = E([0], [g0], [r0])
q0 = Q(0, f0, u0, e0, Dom=any, Avoid=avoid)
# [0.086, 0.086, -1.34 + (-0.38+1.34)/2.0, 0.58+(1.54-0.58)*(3/4)]
# [0.5, -0.3, -3, -2]
# [ 0.12402108 -0.26131208 -1.11743065 -0.82005916]
# [-0.22369904  0.16082207  0.6009526  -1.00018365]
# [-0.22402108, 0.16131208, 0.6, -1.0]
init_X = np.array([-0.22369904, 0.16082207, 0.6009526, -1.00018365])
init_qID = 0

states = [r'\theta_{sw}', r'\theta_{st}',
          r'\dot{\theta}_{sw}', r'\dot{\theta}_{st}']

h = H([q0], Init_X = init_X, Init_qID = init_qID, 
      state_names = states )



