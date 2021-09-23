# -*- coding: utf-8 -*-
"""
Created on Sun October 25 15:45:36 2020
MA 514000 Project 2

@author: Sookyung Lee
"""
import numpy as np
import pandas as pd
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
f = lambda x: 1/(1+6*(x**2))
xi = np.arange(-1, 1.2, 0.1)
yi = f(xi)
coeff = cubicNaturalCoeff(xi, yi)

knots = np.linspace(-1, 1, 5)
xvals = np.array(knots)
yvals = np.array(f(knots))
natural_spline = np.array(list(map(lambda x: cubicSplineEval(x, xi, coeff), knots)))

#Produce table of required values
df = pd.DataFrame({"x": xvals, "f(x)": yvals, "p(x)": natural_spline})
df['f(x) - p(x)'] = df['f(x)'] - df['p(x)']
print(df)

#Graph Natural Cubic Spline Interpolation
plt.plot(knots, natural_spline, label='Natural cubic spline')
plt.plot(xi, yi, '.', label='Data Points')
plt.xlim([-1, 1])
plt.legend(loc='best')
