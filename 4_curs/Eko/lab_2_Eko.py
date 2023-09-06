
from scipy.integrate import odeint

import numpy as np

import plotly

import plotly.graph_objs as go

from plotly.offline import init_notebook_mode

#init_notebook_mode(connected=True)



def model(y, t, k1, k2):
    L, D = y

    dydt = [-k1 * L, k1 * L - k2 * D]

    return dydt


t = np.linspace(0, 12, 200)

k1 = 0.1

k2 = 0.2

D0 = 1

L0 = 20

y0 = [L0, D0]

result_20 = odeint(model, y0, t, args=(k1, k2))

L0 = 30

y0 = [L0, D0]

result_30 = odeint(model, y0, t, args=(k1, k2))

L0 = 40

y0 = [L0, D0]

result_40 = odeint(model, y0, t, args=(k1, k2))

k2 = 0.8

result_k2_08 = odeint(model, y0, t, args=(k1, k2))

def analyt(t, k1, k2, y0):
    return k1 / (k2 - k1) * y0[0] * (np.exp(-k1 * t) - np.exp(-k2 * t)) + y0[1] * np.exp(-k2 * t)

def mv_graphics_stat(graphs_visible, my_title,

                     mv_title_x, mv_title_y,

                     file_name,

                     xaxis_min, xaxis_max,

                     mv_legend_x, mv_legend_anchor):
    fig = go.Figure(data=graphs_visible,

                    layout_xaxis_range=[xaxis_min,

                                        xaxis_max])

    fig.layout.font.family = 'Times'

    fig.layout.font.size = 14

    fig.update_layout(xaxis_title=mv_title_x,

                      yaxis_title=mv_title_y,

                      legend=dict(x=mv_legend_x,

                                  xanchor=

                                  mv_legend_anchor),

                      autosize=False,

                      width=600, height=400)

    fig.update_layout(margin=

                      dict(l=25, r=0, t=0, b=25))

    fig.update_layout(yaxis=

                      dict(autorange="reversed"))

    fig.show()

    fig.write_image("images/" + file_name + ".pdf")