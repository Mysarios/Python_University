import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import rc


# источник мощностью M на высоте H
Stranght = 1.4
high = 100.0

# рассматриваем задачу с постоянными коэффициентами
coef_turb_z = 5.0 #Турбулентность по z (поток)
k0 = 0.5
u = 5.0
coef_turb_y = k0 * u #Турбулентность по х (поток)
y = 0.

uH = u
# Сетка значений
z_max = 300.
x_max = 500.

step_x = 0.25
step_z = 1

x = np.arange(0.001, x_max + step_x, step_x)
z = np.arange(0, z_max + step_z, step_z)

n = len(z) - 1
m = len(x) - 1


Wind_speed = np.ones(n + 1)
coef_turb = np.ones(n + 1)

Wind_speed = Wind_speed * u      #скорость ветра - в примере константа
coef_turb = coef_turb * coef_turb_z   #коэф турбулентности - в примере константа

# прогонка
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
        Y[i] = alfa[i] * Y[i+1] + beta[i]

#Массивчики для прогонки
A = np.zeros(n+1)
C = np.zeros(n+1)
B = np.zeros(n+1)
F = np.zeros(n+1)

#Массивчики для результатов
results = np.zeros((n + 1, m + 1))
results_y0 = np.zeros((n + 1, m + 1))
result_not_explicit = np.zeros((n + 1, m + 1))
result_Krank_Nikolson = np.zeros((n + 1, m + 1))
result_analyt = np.zeros((n + 1, m + 1))

index_first_high =0
for index in range(0,len(z)):
    if z[index] == 100:
        index_first_high = index

# Начальные условия для некоторых методов
results[index_first_high, 0] = Stranght / (uH * step_z)
result_not_explicit[index_first_high, 0] = Stranght / (uH * step_z)
result_Krank_Nikolson[index_first_high, 0] = Stranght / (uH * step_z)


print("Start excplicit result")
for j in range(m):
    for i in range(1, n-1):
        sigma = step_x / (Wind_speed[i] * step_z ** 2)
        k_minus = 0.5 * (coef_turb[i] + coef_turb[i - 1])
        k_plus = 0.5 * (coef_turb[i] + coef_turb[i + 1])
        #c[i, j+1] = (sigma * k_plus * c[i+1, j] + (1- sigma * (k_minus + k_plus)) * c[i, j]
         #           + sigma * k_minus * c[i-1, j])
        results[i, j + 1] = (sigma * (k_plus * results[i + 1, j] + k_minus * results[i - 1, j]) +
                             (1 - sigma * (k_minus + k_plus)) * results[i, j])

    results[0, j + 1] = results[1, j + 1]
    results[n, j] = 0
for j in range(1, m):
    for i in range(n):
        results_y0[i,j] = results[i,j] / (2 * math.sqrt(math.pi * k0 * x[j]))

print("end Result get: ",results_y0)

print("Start not - excplicit result")
for j in range(1, m+1):
    # вычисление коэффциентов системы уравнений
    C[0] = 1.; B[0] = 1.; F[0] = 0.
    C[n] = 1.; A[n] = 0.; F[n] = 0.
    for i in range(1,n):
        sigma = Wind_speed[i] * step_z ** 2 / step_x
        k_minus = 0.5 * (coef_turb[i] + coef_turb[i - 1])
        k_plus = 0.5 * (coef_turb[i] + coef_turb[i + 1])
        A[i] = k_minus
        C[i] = k_plus + k_minus + sigma
        B[i] = k_plus
        F[i] = sigma * result_not_explicit[i, j - 1] #*0.983
    progon(A, C, B, F, result_not_explicit[:, j])

print("end Result get: ",result_not_explicit)

print("Start Krank-Nikolson result")
for j in range(1, m+1):
    # вычисление коэффциентов системы уравнений
    C[0] = 1.; B[0] = 1.; F[0] = 0.
    C[n] = 1.; A[n] = 0.; F[n] = 0.
    for i in range(1,n):
        sigma = Wind_speed[i] * step_z ** 2 / step_x
        k_minus = 0.5 * (coef_turb[i] + coef_turb[i - 1])
        k_plus = 0.5 * (coef_turb[i] + coef_turb[i + 1])
        A[i] = k_minus
        C[i] = k_plus + k_minus + 2 * sigma
        B[i] = k_plus
        F[i] = (2 * sigma * result_Krank_Nikolson[i, j - 1] + k_plus * (result_Krank_Nikolson[i + 1, j - 1] - result_Krank_Nikolson[i, j - 1])
                - k_minus * (result_Krank_Nikolson[i, j - 1] - result_Krank_Nikolson[i - 1, j - 1]))
    progon(A, C, B, F, result_Krank_Nikolson[:, j])

print("end Result get: ", result_Krank_Nikolson)

# сетка для графика
x_array_grid, z_array_grid = np.meshgrid(x, z)

# определение линий уровня и цветов
Line_highs = [0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.016, 0.020, 0.030, 0.040, 0.300]
Color_lines = ['cyan', 'turquoise', 'teal', 'steelblue', 'cornflowerblue', 'mediumslateblue',
       'indigo', 'darkmagenta', 'orchid', 'lightpink', 'hotpink']

temp = np.zeros(m+1)
temp1 = np.zeros(m+1)
temp2 = np.zeros(m+1)
temp3 = np.zeros(m+1)
temp4 = np.zeros(m+1)
temp5 = np.zeros(m+1)
temp6 = np.zeros(m+1)
index_first_high, = np.where(z == 90)

temp[:]  = results[index_first_high, :] / (2 * np.sqrt(math.pi * k0 * x[:])) * np.exp(-y ** 2 / (4 * k0 * x[:]))
temp1[:] = result_not_explicit[index_first_high, :] / (2 * np.sqrt(math.pi * k0 * x[:])) * np.exp(-y ** 2 / (4 * k0 * x[:]))
temp2[:] = (Stranght / (4 * math.pi * x[:] * math.sqrt(coef_turb_y * coef_turb_z))
            * np.exp(-u * y ** 2 / (4 * coef_turb_y * x[:]))
            * (np.exp(-u * (z[index_first_high] + high) ** 2 / (4 * coef_turb_z * x[:])) + np.exp(-u * (z[index_first_high] - high) ** 2 / (4 * coef_turb_z * x[:])))
            )
temp3[:] = result_Krank_Nikolson[index_first_high, :] / (2 * np.sqrt(math.pi * k0 * x[:])) * np.exp(-y ** 2 / (4 * k0 * x[:]))

index_first_high, = np.where(z == 70)
temp4[:] = result_Krank_Nikolson[index_first_high, :] / (2 * np.sqrt(math.pi * k0 * x[:])) * np.exp(-y ** 2 / (4 * k0 * x[:]))
index_first_high, = np.where(z == 80)
temp5[:] = result_Krank_Nikolson[index_first_high, :] / (2 * np.sqrt(math.pi * k0 * x[:])) * np.exp(-y ** 2 / (4 * k0 * x[:]))

index_first_high, = np.where(z == 100)
temp6[:] = result_Krank_Nikolson[index_first_high, :] / (2 * np.sqrt(math.pi * k0 * x[:])) * np.exp(-y ** 2 / (4 * k0 * x[:]))

print('явная    неявная      аналитическое')
for i in range(0, 20):
    print(temp[i], temp1[i], temp2[i], temp3[i])


plt.plot(x, temp, c = "r", label = "явная")
plt.plot(x, temp1, c = "g", label = "неявная")
plt.plot(x, temp2, c = "b", label = "аналитическое")
plt.plot(x, temp3, c = "m", label = "схема Кранка-Николсона")
plt.plot(x, temp4, c = "c", label = "К-Н высота 70")
plt.plot(x, temp5, c = "r", label = "К-Н высота 110")
plt.plot(x, temp6, c = "g", label = "К-Н высота 100")

plt.ylim([0.0, 0.001])
plt.xlim([0.05, 500.])
plt.xlabel('$x$, м')
plt.ylabel('$z$, м')
plt.title('Все графики', loc='left')
plt.legend()
plt.show()
