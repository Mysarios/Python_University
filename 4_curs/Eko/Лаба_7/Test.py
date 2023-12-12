import math
import numpy as np
import matplotlib.pyplot as plt
from numba import njit
import plotly.graph_objects as go
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from plotly.offline import init_notebook_mode
#import kaleido
import plotly
#from tqdm.notebook import tqdm

# Входные параметры
# Расчетная область
# Область определения по x
xmin, xmax, dx = 0., 9000., 30.
# Область определения по y
ymin, ymax, dy = 0., 9000., 30.
# Область определения по z
zmin, zmax, dz = 0., 540., 30.
# Время исследования
t1, t2, dt = 100, 1000, 10

# Сетка по x
x = np.arange(xmin, xmax, dx)
# Сетка по y
y = np.arange(ymin, ymax, dy)
nx = len(x)
ny = len(y)
# Высота источника
hs = 30.
# Координаты источника
xss, yss, zss = 0, 0, hs
# Средние скорости ветра, м/с
#U, V, W = 3., 3., 0.3
U, V, W = 3., 0.0003, 0.003
# Мощность выброса
M = 100.
Mm = M * 1000
# Стратификация
s = 2 # 1, 2, 3
# Коэффициенты турбулентности
i1, i2 = 0.1, 0.1
K1, K2 = 10., 1000.

# Вычисляем K3 в зависимости от стратификации
if s==1:
    i3 = 0.2
    K3 = 50.
if s==2:
    i3 = 0.1
    K3 = 10.
if s==3:
    i3 = 0.05
    K3 = 1.
# Диаметр источника
#ds = 0.8
# Радиус источника
#R0 = ds/2

#Ta, w0, Ts = 298, 10, 358

# Дисперсии атмосферы
sig1 = U * i1
sig2 = V * i2
sig3 = W * i3

# Лагранжевы временные масштабы
T11 = K1 / sig1**2
T12 = K2 / sig2**2
T13 = K3 / sig3**2

p11 = math.exp(-dt/T11)
p12 = math.exp(-dt/T12)
p13 = math.exp(-dt/T13)

# Дисперсии нормального распределения
sigeps1 =  sig1**2 * (1-p11**2)
sigeps2 =  sig2**2 * (1-p12**2)
sigeps3 =  sig3**2 * (1-p13**2)

n = 360
Ccalc = np.zeros((n,n))
CLcalc = np.zeros((n,n,n))


@njit
def calpuff():
    n = 360
    # global
    Cc, Ccc, Ca = np.zeros((n, n)), np.zeros(n), np.zeros((n, n))
    CLagrange = np.zeros((n, n, n))
    xc, C = np.zeros(n), np.zeros((n, n, n))
    eps1, eps2, eps3 = np.zeros(n), np.zeros(n), np.zeros(n)
    Us1, Us2, Us3 = np.zeros(n), np.zeros(n), np.zeros(n)
    U1, U2, U3 = np.zeros(n), np.zeros(n), np.zeros(n)
    sigx1, sigx2, sigx3 = np.zeros(n), np.zeros(n), np.zeros(n)
    xc1, xc2, xc3 = np.zeros(n), np.zeros(n), np.zeros(n)
    #    Cp = np.zeros(n)
    t = 0
    tk = t1
    #    j = 0
    #    tp = 1200
    while True:
        #  for t in tqdm(range(1, t2)):
        t += 1
        # Момент запуска источника
        if t == t1:
            K = 0  # Количество клубов
            #           j +=1
            print('Запуск источника!!!')
            print('t1 = ', t1, ' t = ', t)
            Xx = np.random.randn()
            eps1[K] = np.sqrt(sigeps1) * Xx
            Xx = np.random.randn()
            eps2[K] = np.sqrt(sigeps2) * Xx
            Xx = np.random.randn()
            eps3[K] = np.sqrt(sigeps3) * Xx

            Us1[K], Us2[K], Us3[K] = eps1[K], eps2[K], eps3[K]

            U1[K] = U + Us1[K]
            U2[K] = V + Us2[K]
            U3[K] = W + Us3[K]
            sigx1[K], sigx2[K], sigx3[K] = 0, 0, 0
            xc1[K], xc2[K], xc3[K] = xss, yss, zss

            tk = t1 + dt
            # Работа источника
        if t == tk and t <= t2:
            K += 1
            # print('Работа источника!!!')
            print(' t = ', t)
            Xx = np.random.randn()
            eps1[K] = np.sqrt(sigeps1) * Xx
            Xx = np.random.randn()
            eps2[K] = np.sqrt(sigeps2) * Xx
            Xx = np.random.randn()
            eps3[K] = np.sqrt(sigeps3) * Xx

            Us1[K] = eps1[K]
            Us2[K] = eps2[K]
            Us3[K] = eps3[K]

            U1[K] = U + Us1[K]
            U2[K] = V + Us2[K]
            U3[K] = W + Us3[K]

            sigx1[K] = 0
            sigx2[K] = 0
            sigx3[K] = 0

            xc1[K] = xss
            xc2[K] = yss
            xc3[K] = zss

            for k in range(K):
                Xx = np.random.randn()
                eps1[k] = np.sqrt(sigeps1) * Xx
                Xx = np.random.randn()
                eps2[k] = np.sqrt(sigeps2) * Xx
                Xx = np.random.randn()
                eps3[k] = np.sqrt(sigeps3) * Xx

                Us1[k] = Us1[k] * p11 + eps1[k]
                Us2[k] = Us2[k] * p12 + eps2[k]
                Us3[k] = Us3[k] * p13 + eps3[k]

                U1[k] = U + Us1[k]
                U2[k] = V + Us2[k]
                U3[k] = W + Us3[k]

                td = t - t1 - k * dt

                # print(td)
                sigx1[k] = 2 * K1 * td
                sigx2[k] = 2 * K2 * td
                sigx3[k] = 2 * K3 * td

                xc1[k] += U1[k] * dt
                xc2[k] += U2[k] * dt
                xc3[k] += U3[k] * dt
            # print(xc3[k], hs)

            tk = tk + dt
            #            j += 1

            for i in range(nx):
                for l in range(ny):
                    for nn in range(K):
                        C[i, l, nn] = (2 * Mm * dt / ((2 * math.pi) ** 1.5 * np.sqrt(sigx1[nn] * sigx2[nn] * sigx3[nn]))
                                       * np.exp(-(x[i] - xc1[nn]) ** 2 / (2 * sigx1[nn]))
                                       * np.exp(-(y[l] - xc2[nn]) ** 2 / (2 * sigx2[nn]))
                                       * np.exp(-xc3[nn] ** 2 / (2 * sigx3[nn]))
                                       # * np.exp(-hs**2/(2*sigx3[nn]))
                                       # * np.exp(-xc3[nn]**2/(2*sigx3[nn]))
                                       )

            for i in range(nx):
                for l in range(ny):
                    Cc[i, l] = np.sum(C[i, l, :])
                    #   print(Cc[:,0])
            for i in range(nx):
                Ccc[i] = np.sum(Cc[i, :])
            for i in range(nx):
                for l in range(ny):
                    CLagrange[i, l, K] = Cc[i, l]
        if t == 3700:
            break
    return Cc, CLagrange  # , Ccc\n

Ccalc, CLcalc = calpuff()

xa = np.arange(xmin+0.001, xmax+0.001, dx)
Ca = np.zeros((n,n))
for i in range(1, nx):
    for l in range(ny):
#        Ca[i,l] = (Mm/(2*math.pi*x[i]*np.sqrt(K2*K3))
 #             * np.exp(-U*hs**2/(4*K3*x[i]) - U*y[l]**2/(4*K2*x[i]))
 #               )
        Ca[i,l] = (Mm/(4*math.pi*xa[i]*np.sqrt(K2*K3))
                 * np.exp(-U*y[l]**2/(4*K2*xa[i]))
                * (np.exp(-U*(hs)**2 / (4*K3*xa[i]))+np.exp(-U*(hs)**2 / (4*K3*xa[i])))
      #         * np.exp(-U*y[l]**2/(4*K2*xa[i]))
     #         * np.exp(- U*hs**2/(4*K3*xa[i]))
                )

CLsum = np.zeros((n,n))
for i in range(nx):
    for j in range(ny):
        CLsum[i,j] = CLcalc[i,j,240:359].mean()


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
    fig.write_image(file_name + ".png")
    fig.show()


graphs_visible = [go.Scatter(visible=True, x=x, y=CLsum[:, 0],
                             name='Расчет',
                             line_dash='dot', line_color='#6e0cb3'),
                  go.Scatter(visible=True, x=xa, y=Ca[:, 0],
                             name='Аналитическое решение',
                             line_dash='solid', line_color='#6e0cb3'),

                  ]

my_title = "Приземные концентрации"
mv_title_x = '$$x, м$$'
mv_title_y = '$$c^{\prime}, мг/м^3$$'

file_name = "analit_calc"
xaxis_min = 0.
xaxis_max = 9000.
mv_legend_x = 1.
mv_legend_anchor = "right"

mv_graphics_stat(graphs_visible, my_title,
                 mv_title_x, mv_title_y,
                 file_name,
                 xaxis_min, xaxis_max,
                 mv_legend_x, mv_legend_anchor)

