import math
from numba import njit
import numpy as np
import matplotlib.pyplot as plt
import plotly
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from plotly.offline import init_notebook_mode
init_notebook_mode(connected = True)



#Входные
# Постоянная Кармана
KAPPA = 0.4187
L = 100. # Масштаб
Stranght = 10.
high = 10. / abs(L)

U = 10.
coef_turb = 5.

uH = U
UDIN = 0.5
# Коэффициент для перехода к безразмерной концентрации
C = Stranght / (UDIN * abs(L))

# Сетка
z_step = 0.05
x_step = 0.01

z0 = 0.0002
limitX = 3000.
xmax = limitX / L   # + 0.00001
zmax = 1. + z0


def psi(z, L, z0):
    return 1. / (1. + 4.7 * (z + z0)) if L > 0 else (1 + (1 if L > 0 else -1) * 15 * z) ** 0.25
def iks(z, L):
    return (1 + (1 if L > 0 else -1) * 15 * z) ** 0.25
def f_MO(z, L):
    return (np.log(z * L) + 4.7 * z if L > 0 else
            (np.log(z * abs(L)) - 2. * np.log(0.5 * (1 + iks(z, L))) -
             np.log(0.5 * (1 + iks(z, L))) + 2. * np.arctan(iks(z, L))))
def u_MO(z, L, z0):
    return 1. / KAPPA * (f_MO(z, L) - f_MO(z0, L))
def k_MO(z, L, z0):
    return KAPPA * (z + z0) * psi(z, L, z0)
@njit
def progon(A, C, B, F, Y):
    alfa = np.zeros(n + 1)
    beta = np.zeros(n + 1)

    # прямой ход метода прогонки
    # вычисление коэффициентов прогонки
    alfa[0] = B[0] / C[0]
    beta[0] = F[0] / C[0]
    for i in range(1, n + 1):
        alfa[i] = B[i] / (C[i] - alfa[i - 1] * A[i])
        beta[i] = (A[i] * beta[i - 1] + F[i]) / (C[i] - alfa[i - 1] * A[i])

    # обратный ход метода прогонки
    Y[i] = beta[n]
    for i in range(n - 1, -1, -1):
        Y[i] = alfa[i] * Y[i + 1] + beta[i]
@njit
def explicit_scheme():
    # вычисление по явной схеме
    for j in range(m):
        for i in range(1, n - 1):
            sigma = x_step / (Wind_speed_main[i] * z_step ** 2)
            k_minus = 0.5 * (Coef_turb_main[i] + Coef_turb_main[i - 1])
            k_plus = 0.5 * (Coef_turb_main[i] + Coef_turb_main[i + 1])
            # c[i, j+1] = (sigma * k_plus * c[i+1, j] + (1- sigma * (k_minus + k_plus)) * c[i, j]
            #           + sigma * k_minus * c[i-1, j])
            result[i, j + 1] = (sigma * (k_plus * result[i + 1, j] + k_minus * result[i - 1, j]) +
                           (1 - sigma * (k_minus + k_plus)) * result[i, j])

        result[0, j + 1] = result[1, j + 1]
        result[n, j] = 0

    for j in range(1, m):
        for i in range(n):
            cpoint[i, j] = result[i, j] / (2 * math.sqrt(math.pi * k0 * x[j]))
    # вычисление по явной схеме выполнено
@njit
def implicit_scheme(cn):
    # global cn
    # определение массивов для коэффициентов системы
    # для метода прогонки в чисто неявной схеме
    A = np.zeros(n + 1)
    C = np.zeros(n + 1)
    B = np.zeros(n + 1)
    F = np.zeros(n + 1)
    # вычисление по чисто неявной схеме
    for j in range(1, m + 1):
        # вычисление коэффциентов системы уравнений
        C[0], B[0], F[0] = 1., 1., 0.
        C[n], A[n], F[n] = 1., 0., 0.
        for i in range(1, n):
            sigma = Wind_speed_main[i] * z_step ** 2 / x_step
            k_minus = 0.5 * (Coef_turb_main[i] + Coef_turb_main[i - 1])
            k_plus = 0.5 * (Coef_turb_main[i] + Coef_turb_main[i + 1])
            A[i] = k_minus
            C[i] = k_plus + k_minus + sigma
            B[i] = k_plus
            F[i] = sigma * cn[i, j - 1]  # *0.983
        progon(A, C, B, F, cn[:, j])
    # вычисление по неявной схеме выполнено
@njit
def Krank_Nikolson():
    # вычисление по схеме Кранка-Николсона
    for j in range(1, m + 1):
        # вычисление коэффциентов системы уравнений
        C[0] = 1.;
        B[0] = 1.;
        F[0] = 0.
        C[n] = 1.;
        A[n] = 0.;
        F[n] = 0.
        for i in range(1, n):
            sigma = Wind_speed_main[i] * z_step ** 2 / x_step
            k_minus = 0.5 * (Coef_turb_main[i] + Coef_turb_main[i - 1])
            k_plus = 0.5 * (Coef_turb_main[i] + Coef_turb_main[i + 1])
            A[i] = k_minus
            C[i] = k_plus + k_minus + 2 * sigma
            B[i] = k_plus
            F[i] = (2 * sigma * ck[i, j - 1] + k_plus * (ck[i + 1, j - 1] - ck[i, j - 1])
                    - k_minus * (ck[i, j - 1] - ck[i - 1, j - 1]))

        progon(A, C, B, F, ck[:, j])
def mv_graphics_stat(graphs_visible, my_title,mv_title_x, mv_title_y,file_name,xaxis_min, xaxis_max,mv_legend_x, mv_legend_anchor):
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

    fig.show()


x = np.arange(0., xmax + x_step, x_step)
z = np.arange(z0, zmax + z_step, z_step)

n = len(z) - 1
m = len(x) - 1
result = np.zeros((n + 1, m + 1))

indexH=0
for i in range(0,len(z)):
    if z[i-1] <= high and high <= z[i]:
        indexH = i
        break

Coef_turb_main = np.zeros(n)
Wind_speed_main = np.zeros(n)
Coef_turb_main = k_MO(z, L, z0)
Wind_speed_main = u_MO(z, L, z0)

result[indexH, 0] = 1. / (Wind_speed_main[indexH] * z_step)  #Начальные условия
implicit_scheme(result)

# вычисление по схеме Кранка-Николсона выполнено
graphs_visible = [go.Scatter(visible = True, x=x*L, y=result[1, :] * C, name='{:,.1f} м'.format(z[1]*L).replace(",", " "),
                             line_dash = 'dot', line_color = '#ffa306'),
                    go.Scatter(visible = True, x=x*L, y=result[4, :] * C,
                               name='{:,.1f} м'.format(z[4]*L).replace(",", " "),
                               line_dash = 'dot', line_color = '#ffa306'),
                    go.Scatter(visible = True, x=x*L, y=result[8, :] * C,
                               name='{:,.0f} м'.format(z[8]*L).replace(",", " "),
                               line_dash = 'solid', line_color = '#ffa306'),
                    go.Scatter(visible = True, x=x*L, y=result[12, :] * C,
                               name='{:,.0f} м'.format(z[12]*L).replace(",", " "),
                               line_dash = 'dot', line_color = '#2a11ff'),
                 ]

my_title = "Приземные концентрации на разных высотах"
mv_title_x = "$$t, ч$$"
mv_title_y = "$$c,  ^\circ С$$"
file_name = "diary_tau"
xaxis_min = 0.
xaxis_max = 3000.
mv_legend_x = 1.
mv_legend_anchor = "right"

mv_graphics_stat(graphs_visible, my_title,
                 mv_title_x, mv_title_y,
                 file_name,
                 xaxis_min, xaxis_max,
                 mv_legend_x, mv_legend_anchor)

graphs_visible = [go.Scatter(visible=True, x=Wind_speed_main * UDIN, y=z,
                             name='{:,.1f} м'.format(z[4]).replace(",", " "),
                             line_dash = 'dot', line_color = '#ffa306')
                 ]

my_title = "Приземные концентрации на разных высотах"
mv_title_x = "$$u, м/с$$"
mv_title_y = "$$z,  м$$"
file_name = "diary_tau"
xaxis_min = 0.
xaxis_max = 33.
mv_legend_x = 1.
mv_legend_anchor = "right"

mv_graphics_stat(graphs_visible, my_title,
                 mv_title_x, mv_title_y,
                 file_name,
                 xaxis_min, xaxis_max,
                 mv_legend_x, mv_legend_anchor)

graphs_visible = [go.Scatter(visible=True, x=Coef_turb_main, y=z,
                             name='{:,.1f} м'.format(z[4]).replace(",", " "),
                             line_dash = 'dot', line_color = '#ffa306')
                 ]


my_title = "Приземные концентрации на разных высотах"
mv_title_x = "$$k, $$"
mv_title_y = "$$z, м$$"
file_name = "diary_tau"
xaxis_min = 0.
xaxis_max = 0.1
mv_legend_x = 1.
mv_legend_anchor = "right"

mv_graphics_stat(graphs_visible, my_title,
                 mv_title_x, mv_title_y,
                 file_name,
                 xaxis_min, xaxis_max,
                 mv_legend_x, mv_legend_anchor)