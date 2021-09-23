# -*- coding: utf-8 -*-
"""
Created on Sun October 25 12:22:41 2020
MA 514000 Project 2

@author: Sookyung Lee
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def divided_difference(x_vals, y_vals):
    datapts = x_vals.size
    c_vals = np.zeros(datapts)

    for i in range(datapts):
        c_vals[i] = y_vals[i]

    for i in range(1, datapts):
        for j in np.flip(np.arange(i, datapts)):
            d_vals = x_vals[j] - x_vals[j - i]
            u_vals = c_vals[j] - c_vals[j - 1]

            c_vals[j] = u_vals/d_vals

    return c_vals


def newton(x_coor, y_coor, z_poly):
    datapts = x_coor.size
    diff = divided_difference(x_coor, y_coor)

    u_val = diff[0]
    c_vals = 1.0

    for i in range(datapts - 1):
        c_vals *= (z_poly- x_coor[i])
        u_val += diff[i + 1]*c_vals
    return u_val



def chebyshev(i):
    return np.cos((i-1)*np.pi/20)



#Running the numerical experiment
f = lambda x: 1/(1+6*(x**2))
i = np.linspace(1, 21, 21)
xi = chebyshev(i[:])
yi = f(xi)


nodes = np.linspace(-1, 1, 41)
fvalue = f(nodes)

interp = newton(xi, yi, nodes)
xvals = np.array(nodes)
yvals = np.array(f(nodes))

interpolation = np.array(interp)
df = pd.DataFrame({"x": xvals, "f(x)": yvals, "p(x)": interpolation})
df['f(x) - p(x)'] = df['f(x)'] - df['p(x)']
print(df)

#Graph Interpolation
plt.plot(nodes, interpolation, label='Interpolation with Chebyshev')
plt.plot(xi, yi, '.', label='Data Points')
plt.xlim([-1, 1])
plt.legend(loc='best')
