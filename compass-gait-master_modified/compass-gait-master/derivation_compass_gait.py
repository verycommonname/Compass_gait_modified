#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Rafael Figueroa
Compass-Gait Dynamics derivation from Gaswani
Gaswani equations confirmed by first principles derivation
using Lagrangian
"""

from __future__ import division
import numpy as np
from sympy import *
import shelve
init_printing()

m, l, g, M, a, b = symbols('m, l, g, M, a, b')
hst, hsw, hdst, hdsw = symbols('hst, hsw, hdst, hdsw')
u = symbols('u')

# Dynamics on Standard Manipulator Form
H00 = m*b**2
H01 = -m*l*b*cos(hst-hsw)
H10 = -m*l*b*cos(hst-hsw)
H11 = (M+m)*l**2+m*a**2

C00 = 0
C01 = m*l*b*sin(hst-hsw)*hdst
C10 = -m*l*b*sin(hst-hsw)*hdsw
C11 = 0

G00 = m*b*g*sin(hsw)
G10 = -(M*l+m*a+m*l)*g*sin(hst)

B00 = 1.0
B10 = -1.0

H = Matrix([[H00, H01],[H10, H11]])
C = Matrix([[C00, C01],[C10, C11]])
G = Matrix(2, 1, [G00, G10])
B = Matrix(2, 1, [B00, B10])

qd = Matrix(2, 1, [hdsw, hdst])

# Transition at impact calculation
alpha = (1/2)*(hsw-hst)

Q_before00 = -m*a*b
Q_before01 = -m*a*b+(M*l**2+2*m*a*l)*cos(2*alpha)
Q_before10 = 0
Q_before11 = -m*a*b

Q_after00 = m*b*(b-l*cos(2*alpha))
Q_after01 = m*l*(l-b*cos(2*alpha)) + m*a**2 + M*l**2
Q_after10 = m*b**2
Q_after11 = -m*b*l*cos(2*alpha)

Q_before = Matrix([[Q_before00, Q_before01],[Q_before10, Q_before11]])
Q_after = Matrix([[Q_after00, Q_after01],[Q_after10, Q_after11]])

P = simplify(Q_after.inv() * Q_before)

# Pretty printing

print ('H')
pprint(H)
print ('C')
pprint(C)
print ('G')
pprint(G)
print ('B')
pprint(B)
print('qd')
pprint(qd)
print ('C*qd')
pprint(C*qd)
print ('B*u')
pprint(B*u)

qdd = simplify(H.inv()*(-C*qd -G + B*u))
print('qdd')
pprint(qdd)

qdd0 = qdd[0]
qdd1 = qdd[1]
print('qdd 0 and 1')
pprint(qdd0)
pprint(qdd1)

# Expressions for use in the dynamics functions at simulation
print('Expression for swing leg:')
print(qdd0)
print('Expression for standing leg:')
print(qdd1)
print('Expression for P')
print(P)

# Pointcare Map Derivation

W = Matrix([[0, 1, 0, 0],
            [1, 0, 0, 0],
            [0, 0, P[0, 0], P[0, 1]],
            [0, 0, P[1, 0], P[1, 1]]])

print ('P and W')
pprint(P)
pprint(W)

A = Matrix([[0, 0, 1, 0],
            [0, 0, 0, 1],
            [-23.91, 25.73, 0, 0],
            [-4.54, 15.44, 0, 0]])

T = symbols('T') 
stm = mpmath.expm(A*T)
print ('State Transition Matrix')
stmSym = Matrix(stm)
pprint(stmSym)

D = dot(W, stmSym)
pprint(D)






