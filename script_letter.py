# -*- coding: utf-8 -*-
"""
Created on Sun October 25 12:15:06 2020
MA 514000 Project 2

@author: Sookyung
"""
import numpy as np
import matplotlib.pyplot as plt

def cubicNaturalCoeff(x_val, y_val):
    datapts = x_val.size
    end_idx = datapts - 1

    a_coeff = np.zeros(datapts)
    b_coeff = np.zeros(end_idx)
    c_coeff = np.zeros(datapts)
    d_coeff = np.zeros(end_idx)

    for i in range(datapts):
        a_coeff[i] = y_val[i]

    h_val = np.zeros(end_idx)

    for i in range(end_idx):
        h_val[i] = x_val[i + 1] - x_val[i]

    u_val = np.zeros(end_idx)
    u_val[0] = 0

    for i in range(1, end_idx):
        u_val[i] = 3*(a_coeff[i + 1] - a_coeff[i])/h_val[i]-3*(a_coeff[i] - a_coeff[i - 1])/h_val[i - 1]
    s_val = np.zeros(datapts)
    z_val = np.zeros(datapts)
    t_val = np.zeros(end_idx)
    s_val[0] = 1
    z_val[0] = 0
    t_val[0] = 0

    for i in range(1, end_idx):
        s_val[i] = 2*(x_val[i + 1] - x_val[i - 1]) - h_val[i - 1] * t_val[i - 1]
        t_val[i] = h_val[i]/s_val[i]
        z_val[i] = (u_val[i] - h_val[i - 1]*z_val[i - 1])/s_val[i]

    s_val[datapts-1] = 1
    z_val[datapts-1] = 0
    c_coeff[datapts-1] = 0

    for i in np.flip(np.arange(end_idx)):
        c_coeff[i] = z_val[i]-t_val[i]*c_coeff[i+1]
        b_coeff[i] = (a_coeff[i+1]-a_coeff[i])/h_val[i]-h_val[i]*(c_coeff[i+1]+2*c_coeff[i])/3
        d_coeff[i] = (c_coeff[i+1]-c_coeff[i])/(3*h_val[i])
    return a_coeff, b_coeff, c_coeff, d_coeff


def cubicSplineEval(t, x_val, coeff):
    datapts = x_val.size
    a_coeff = coeff[0]
    b_coeff = coeff[1]
    c_coeff = coeff[2]
    d_coeff = coeff[3]

    #check the bounds
    if t < x_val[0] or t > x_val[datapts -1]:
        return

    end_idx = datapts - 1
    subint = 0
    for i in range(end_idx):
        if t <= x_val[i + 1]:
            break
        else:
            subint += 1
            
    eval = a_coeff[subint]+b_coeff[subint]*(t - x_val[subint])+c_coeff[subint]*(t-x_val[subint])**2+d_coeff[subint]*(t-x_val[subint])**3

    return eval


#Run the numerical experiment
t = np.linspace(1, 11, 11)
xi = np.array([3.1, 3.05, 2.20, 2.20, 3.6 , 6.2, 6.8, 6.2, 4.4, 4, 4.4])
yi = np.array([6.8, 4.5, 0.9, 2.7, 1.7, 2.2, 5.3, 8.6, 10, 5.3, 0.5])
sx_coeff = cubicNaturalCoeff(t, xi)
sy_coeff = cubicNaturalCoeff(t, yi)

points = np.linspace(0, 6000, 1000000)

sx_naturalspline = np.array(list(map(lambda x: cubicSplineEval(x, t, sx_coeff), points)))
sy_naturalspline = np.array(list(map(lambda x: cubicSplineEval(x, t, sy_coeff), points)))


plt.plot(sx_naturalspline, sy_naturalspline, label='Natural cubic spline')
plt.plot(xi, yi, 'o', label='Data Points')
plt.xlim([0, 10])
plt.legend(loc='best')


