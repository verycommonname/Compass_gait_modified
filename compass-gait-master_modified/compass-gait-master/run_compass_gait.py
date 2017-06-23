#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Rafael Figueroa
Run Compass Gait Simulation
"""

from __future__ import division
from robot_vars import *
import numpy as np
from hasimpy import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import compass_gait_model
from compass_gait_model import h

sin = np.sin
cos = np.cos
tan = np.tan
pi = np.pi

# print all the array
np.set_printoptions(threshold=np.inf)


def psiZ(z):
    zPsi = ((z))/sum(((z)))
    return zPsi


def show_plots():
    simResult.simPlot()
    simResult.phasePlot([0, 2])
    simResult.phasePlot([1, 3])
    plt.show()
    input('\n Press ENTER to close plots') # raw_input()


def point0(hst, hsw):
    """Standing leg foot position and origin"""
    p0x = 0
    p0x = 0
    return np.array([p0x, p0x])


def point3(hst, hsw):
    """Hip position"""
    p3x = -l*sin(hst)
    p3y =  l*cos(hst)
    return np.array([p3x, p3y])


def point4(hst, hsw):
    """Swing leg foot position"""
    p4x = l*(sin(hsw)-sin(hst))
    p4y = l*(cos(hst)-cos(hsw))
    return np.array([p4x, p4y])


def update_line(frame_id):
    global leg1, leg2

    t = np.arange(0.0, 2.0, 0.01)
    s = np.sin(2 * np.pi * t)
    """Updates leg positions on existing animation plot"""
    sim_X = simPath[frame_id]

    hsw = sim_X[0]
    hst = sim_X[1]

    p_st = point0(hst, hsw)
    p_hp = point3(hst, hsw)
    p_sw = point4(hst, hsw)
    
    # First for standing leg
    leg1.set_data([[p_st[0], p_hp[0]], [p_st[1], p_hp[1]]])
    leg2.set_data([[p_hp[0], p_sw[0]], [p_hp[1], p_sw[1]]])

    return [leg1, leg2] 


def animation_init():
    """Clears the lines from the screen"""
    global leg1, leg2
    leg1.set_data([], [])
    leg2.set_data([], [])
    return [leg1, leg2]


def animate_walking():
    plt.xlim(-2, 2)
    plt.ylim(-2, 2)
    plt.xlabel('x')
    plt.title('leg test')
    # interval in ms, Ts in seconds
    # running 100 times faster
    legs_ani = animation.FuncAnimation(fig1, update_line, len(simPath),
                                       init_func = animation_init,
                                       interval = 0.0, 
                                       blit = True,
                                       repeat = False)

    plt.show()

    input('\n Press ENTER to close plots') # raw_input()


u0 = compass_gait_model.u0
init_X = compass_gait_model.init_X

simResult = h.sim(qID = 0, X = init_X, u = u0, t0 = 0, tlim = 3, debug_flag = False, Ts = 1e-3)

simPath = simResult.path
show_plots()
fig1 = plt.figure()
leg1, = plt.plot([], [], 'r-', lw=3)
leg2, = plt.plot([], [], 'b-', lw=3)
ground, = plt.plot([- 2*l*cos(gamma),
                    2*l*cos(gamma)], 
                    [2*l*sin(gamma), 
                    -2*l*sin(gamma)], 'k-', lw=5)



animation_init()
animate_walking()

