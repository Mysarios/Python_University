# Модель Лотки-Вальтерры

# Линии уровня и поверхность

import numpy as np

import plotly.graph_objects as go

import matplotlib.pyplot as plt
import matplotlib.cm as cm

from plotly.offline import init_notebook_mode

def f(v, a, b3_b1):
    return  b3_b1 * v + a - np.log(a * v**(b3_b1))


V = np.arange(0.01, 5, 0.05)
A = np.arange(0.01, 5, 0.05)


# Creating 2-D grid of features

[X, Y] = np.meshgrid(V,A)


b3_b1 = 1
Z = f(X, Y, b3_b1)
fig,ax = plt.subplots()
plt.contour(X, Y, Z, colors='green');
levels = np.arange(-1.2, 1.6, 0.2)
CS = ax.contour(Z, levels, origin='lower', cmap='flag', extend='both',
                linewidths=2, extent=(-3, 3, -2, 2))
CS.collections[6].set_linewidth(4)

im = ax.imshow(Z, interpolation='bilinear', origin='lower',
               cmap=cm.gray, extent=(0, 5, 0, 5))

CB = fig.colorbar(CS, shrink=0.8)
l, b, w, h = ax.get_position().bounds
ll, bb, ww, hh = CB.ax.get_position().bounds
CB.ax.set_position([ll, b + 0.1*h, ww, h*0.8])
CBI = fig.colorbar(im, orientation='horizontal', shrink=0.8)

plt.show()