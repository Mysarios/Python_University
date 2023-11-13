from sympy import *
import sympy as sym
import numpy as np
import math as m
import matplotlib.pyplot as plt
import scipy.integrate as spi

# Data for programm
A_Numeric = 1
Change = 3
Only_Value_in_Break_point = 1

# Data for algorithm
Q_start = 1.34 / 10
if Only_Value_in_Break_point == 1:
    Q_start = 0

eps = 0.00000005
Max_q_T = 4
N_x = 2
N_y = 2
N = N_x * N_y
Q_max = 245 / (10 ** 6)  # 1 кг максимум

# graph points
Size = 50

# Data about object main
H_coef = 1
h = 0.09 * H_coef
A_lenght_x = 60 * h
B_lenght_y = 60 * h
R_1_x = 225 * h
R_2_y = 225 * h
A = 1
B = 1
K_x = 1 / R_1_x
K_y = 1 / R_2_y
z = 0

# Data material
V = 0.2 * 0.2 * 0.00022  # m^3
P1 = 1500
P2 = 7800
P3 = 1190
Square = 20 * 20
k = 5 / 6
mass = 0

if Change == 0:
    mass = 0
if Change == 1:
    mass = P1 * V
if Change == 2:
    mass = P2 * V
if Change == 3:
    mass = P3 * V
Q_my = (mass / Square) / (10 ** 6)

# Symbols
Xx = Symbol('x')
Yy = Symbol('y')

# Static_DATA for variants:

# Ort:
# T300/976
E1_T300 = 1.4 * (10 ** 5)
E2_T300 = 0.97 * (10 ** 4)
nu_12_T300 = 0.29
nu_21_T300 = (E2_T300 * nu_12_T300) / E1_T300
G_12_T300 = 0.55 * (10 ** 4)
G_13_T300 = G_12_T300
G_23_T300 = 0.33 * (10 ** 4)
F1_min_T300 = -1599
F1_max_T300 = 1517
F2_min_T300 = -253
F2_max_T300 = 46
F12_max_T300 = 41.4
Density_T300 = 1500
Sigma_t_T300 = 100

# Izo:
# Org_Glass_1:
E1_OG = 2.1 * (10 ** 5)
E2_OG = E1_OG
nu_OG = 0.30
G_12_OG = 0.012 * (10 ** 5)
G_13_OG = G_12_T300
G_23_OG = G_12_T300
Sigma_t_OG = 75
Density_OG = 1190

# Still:
E1_Still = 2.1 * (10 ** 5)
E2_Still = E1_Still
nu_Still = 0.3
#G_12_Still = 0.807 * (10 ** 5)
G_12_Still = E1_Still/(2*(1+nu_Still))
G_13_Still = G_12_Still
G_23_Still = G_12_Still
Sigma_t_Still = 265
Density_Still = 7800

# Settings to integral
Start_integral = 0


def Get_W_Plane(x_val, y_val, function, Values, type):
    W_result = 0

    # Hard
    if type == 1:
        for j in range(1, N_x + 1):
            for i in range(1, N_y + 1):
                W_result += Values[(j - 1) * N_x + i - 1] * sin(2 * i * x_val * m.pi / A_lenght_x) * sin(
                    (2 * j - 1) * y_val * m.pi / B_lenght_y)
    if type == 2:
        for j in range(1, N_x + 1):
            for i in range(1, N_y + 1):
                W_result += Values[(j - 1) * N_x + i - 1] * sin((2 * i - 1) * x_val * m.pi / A_lenght_x) * sin(
                    2 * j * y_val * m.pi / B_lenght_y)
    if type == 3:
        for j in range(1, N_x + 1):
            for i in range(1, N_y + 1):
                W_result += Values[(j - 1) * N_x + i - 1] * sin((2 * i - 1) * x_val * m.pi / A_lenght_x) * sin(
                    (2 * j - 1) * y_val * m.pi / B_lenght_y)

    return W_result


def Draw_3d_W(Function, Values_Result, type_f):
    x_array = []
    y_array = []
    z_array = [0] * Size
    Max_Value = 0
    # print(Values_Result)

    # print()
    # print(Function)

    for i in range(0, Size):
        z_array[i] = [0] * Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0, Size):
        x_array.append(i * step_x)
        y_array.append(i * step_y)

    for j in range(0, Size):
        for i in range(0, Size):
            z_array[i][j] = (Get_W_Plane(x_array[i], y_array[j], Function, Values_Result, type_f))
            if abs(z_array[i][j]) > abs(Max_Value):
                Max_Value = z_array[i][j]

    z_array = np.array(z_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x [м]')
    ax.set_ylabel('y [м]')
    ax.set_zlabel(Function + ' (x,y) [м]')
    # ax.set_title(Function + ' (x,y)')
    plt.show()
    return Max_Value


def Draw_3d_Q(Function, Values_Result, type_f, q, Q_function):
    x_array = []
    y_array = []
    z_array = [0] * Size
    Max_Value = 0
    # print(Values_Result)

    # print()
    # print(Function)

    for i in range(0, Size):
        z_array[i] = [0] * Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0, Size):
        x_array.append(i * step_x)
        y_array.append(i * step_y)

    for j in range(0, Size):
        for i in range(0, Size):
            z_array[i][j] = Q_function.subs([(Xx, x_array[i]), (Yy, y_array[j])])

    z_array = np.array(z_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x [м]')
    ax.set_ylabel('y [м]')
    ax.set_zlabel(Function + ' (x,y) [МПа]')
    # ax.set_title(Function + ' (x,y)')
    plt.show()


def Draw_3d_Sigmas(Function, Values_Result, Type_Sigmas, U_function, V_function, W_Function, z_val, W_val, U_val,
                   V_val):
    x_array = []
    y_array = []
    z_array = [0] * Size
    Max_Value = 0

    for i in range(0, Size):
        z_array[i] = [0] * Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0, Size):
        x_array.append(i * step_x)
        y_array.append(i * step_y)

    Sigma_x_Function = 5 * Xx
    Sigma_y_Function = 5 * Xx
    Tay_xy_Function = 5 * Xx

    if Change == 1:
        if Type_Sigmas == 1:
            Sigma_x_Function = (
                Get_Sigma_x_Orto(U_function, V_function, W_Function, E1_T300, nu_12_T300, nu_21_T300, z_val))
        if Type_Sigmas == 2:
            Sigma_y_Function = (
                Get_Sigma_y_Orto(U_function, V_function, W_Function, E2_T300, nu_12_T300, nu_21_T300, z_val))
        if Type_Sigmas == 3:
            Tay_xy_Function = (Get_Sigma_tay_Orto(U_function, V_function, W_Function, G_12_T300, z_val))
    if Change == 2:
        if Type_Sigmas == 1:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
        if Type_Sigmas == 2:
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
        if Type_Sigmas == 3:
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
    if Change == 3:
        if Type_Sigmas == 1:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        if Type_Sigmas == 2:
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        if Type_Sigmas == 3:
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))

    for i in range(N + 1):
        Sigma_x_Function = Sigma_x_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        Sigma_y_Function = Sigma_y_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        Tay_xy_Function = Tay_xy_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])

    # print(Sigma_x_Function)
    # print(Sigma_y_Function)
    # print(Tay_xy_Function)

    for j in range(0, Size):
        for i in range(0, Size):
            if Type_Sigmas == 1:
                z_array[i][j] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if Type_Sigmas == 2:
                z_array[i][j] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if Type_Sigmas == 3:
                z_array[i][j] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if abs(z_array[i][j]) > abs(Max_Value):
                Max_Value = z_array[i][j]

    z_array = np.array(z_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x [м]')
    ax.set_ylabel('y [м]')
    ax.set_zlabel(Function + ' (x,y) [МПа]')
    # ax.set_title(Function + ' (x,y)')
    plt.show()
    return Max_Value


def Draw_3d_Sigmas_main(Function, Values_Result, Type_Sigmas, U_function, V_function, W_Function,PsiX_function, PsiY_function, z_val, W_val, U_val,V_val, PsiX_val, PsiY_val):
    x_array = []
    y_array = []
    z_array = [0] * Size
    Max_Value = 0

    for i in range(0, Size):
        z_array[i] = [0] * Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0, Size):
        x_array.append(i * step_x)
        y_array.append(i * step_y)

    Sigma_x_Function = 5 * Xx
    Sigma_y_Function = 5 * Xx
    Tay_xy_Function = 5 * Xx
    Max_Sigmas_values = [0]*3



    if Change == 1:
        Sigma_x_Function = (Get_Sigma_x_Orto(U_function,V_function,W_Function,E1_T300,nu_12_T300, nu_21_T300,z_val))
        Sigma_y_Function = (Get_Sigma_y_Orto(U_function,V_function,W_Function,E2_T300,nu_12_T300, nu_21_T300,z_val))
        Tay_xy_Function = (Get_Sigma_tay_Orto(U_function,V_function,W_Function,G_12_T300,z_val))
    if Change == 2:
        Sigma_x_Function = Get_Sigma_x_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_OG, nu_OG,
                                           z_val)
        Sigma_y_Function = Get_Sigma_y_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_OG, nu_OG,
                                           z_val)
        Tay_xy_Function = Get_Sigma_tay_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_OG, nu_OG,
                                           z_val)
    if Change == 3:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_Still, nu_Still, z_val))
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_Still, nu_Still, z_val))
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_Still, nu_Still, z_val))


    for i in range(N+1):
        Sigma_x_Function = Sigma_x_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1]),
             ('PsiX' + str(i), PsiX_val[i - 1]), ('PsiY' + str(i), PsiY_val[i - 1])])
        Sigma_y_Function = Sigma_y_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1]),
             ('PsiX' + str(i), PsiX_val[i - 1]), ('PsiY' + str(i), PsiY_val[i - 1])])
        Tay_xy_Function = Tay_xy_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1]),
             ('PsiX' + str(i), PsiX_val[i - 1]), ('PsiY' + str(i), PsiY_val[i - 1])])

    #print(Sigma_x_Function)
    #print(Sigma_y_Function)
    #print(Tay_xy_Function)


    for j in range(0, Size):
        for i in range(0, Size):
                Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                if Change == 1:
                    z_array[i][j] =((Max_Sigmas_values[0] ** 2) + (Max_Sigmas_values[1] ** 2) - Max_Sigmas_values[0] * Max_Sigmas_values[1] + 3 * (Max_Sigmas_values[2] ** 2)) ** (1 / 2)
                if Change == 2:
                    z_array[i][j] =((Max_Sigmas_values[0] ** 2) + (Max_Sigmas_values[1] ** 2) - Max_Sigmas_values[0] * Max_Sigmas_values[1] + 3 * (Max_Sigmas_values[2] ** 2)) ** (1 / 2)
                if Change == 3:
                    z_array[i][j] =((Max_Sigmas_values[0] ** 2) + (Max_Sigmas_values[1] ** 2) - Max_Sigmas_values[0] * Max_Sigmas_values[1] + 3 * (Max_Sigmas_values[2] ** 2)) ** (1 / 2)

                if abs(z_array[i][j]) > abs(Max_Value):
                    Max_Value = z_array[i][j]



    #z_array = np.array(z_array)
    #x_array = np.array(x_array)
    #y_array = np.array(y_array)
    #X, Y = np.meshgrid(x_array, y_array)

    #fig = plt.figure()
    #ax = plt.axes(projection='3d')
    #ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    #ax.set_xlabel('x [м]')
    #ax.set_ylabel('y [м]')
    #ax.set_zlabel(Function + ' (x,y) [МПа]')
    # ax.set_title(Function + ' (x,y)')
    #plt.show()
    print("Max value Miz= ",Max_Value)
    Sigma_Tay_krit = 150
    if Change ==3:
        Sigma_Tay_krit = 265
    return Max_Value/Sigma_Tay_krit


def Draw_3d_Sigmas_main_Orto(Function, Values_Result, Type_Sigmas, U_function, V_function, W_Function, z_val, W_val,
                             U_val, V_val):
    x_array = []
    y_array = []
    z_array = [0] * Size
    Max_Value = 0

    for i in range(0, Size):
        z_array[i] = [0] * Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0, Size):
        x_array.append(i * step_x)
        y_array.append(i * step_y)

    Sigma_x_Function = 5 * Xx
    Sigma_y_Function = 5 * Xx
    Tay_xy_Function = 5 * Xx
    Max_Sigmas_values = [0] * 3

    if Change == 1:
        Sigma_x_Function = (
            Get_Sigma_x_Orto(U_function, V_function, W_Function, E1_T300, nu_12_T300, nu_21_T300, z_val))
        Sigma_y_Function = (
            Get_Sigma_y_Orto(U_function, V_function, W_Function, E2_T300, nu_12_T300, nu_21_T300, z_val))
        Tay_xy_Function = (Get_Sigma_tay_Orto(U_function, V_function, W_Function, G_12_T300, z_val))
    if Change == 2:
        Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
        Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
        Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
    if Change == 3:
        Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))

    for i in range(N + 1):
        Sigma_x_Function = Sigma_x_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        Sigma_y_Function = Sigma_y_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        Tay_xy_Function = Tay_xy_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])

    # print(Sigma_x_Function)
    # print(Sigma_y_Function)
    # print(Tay_xy_Function)
    Buf_1 = (1 / F1_max_T300) + (1 / F1_min_T300)
    Buf_2 = (1 / F2_max_T300) + (1 / F2_min_T300)
    Buf_3 = 1 / (F1_max_T300 * F1_min_T300)
    Buf_4 = 1 / (F2_max_T300 * F2_min_T300)
    Buf_6 = 1 / (F12_max_T300 * F12_max_T300)

    First_val = (1 / 2) * (Buf_4 - Buf_3)
    Second_Val = (1 / 2) * (Buf_3 - Buf_4)
    Third_Val = (1 / 2) * (Buf_3 + Buf_4)

    for j in range(0, Size):
        for i in range(0, Size):
            Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            z_array[i][j] = First_val * (Max_Sigmas_values[0] ** 2) + Second_Val * (Max_Sigmas_values[1] ** 2) \
                            + Buf_1 * Max_Sigmas_values[0] + Buf_2 * Max_Sigmas_values[1] \
                            - Third_Val * ((Max_Sigmas_values[0] - Max_Sigmas_values[2]) ** 2) - Buf_6 * (
                                        Max_Sigmas_values[2] ** 2)
            if abs(z_array[i][j]) > abs(Max_Value):
                Max_Value = z_array[i][j]

    if Only_Value_in_Break_point:
        z_array = np.array(z_array)
        x_array = np.array(x_array)
        y_array = np.array(y_array)
        X, Y = np.meshgrid(x_array, y_array)

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
        ax.set_xlabel('x [м]')
        ax.set_ylabel('y [м]')
        ax.set_zlabel(Function + ' (x,y)')
        # ax.set_title(Function + ' (x,y)')
        plt.show()
    return Max_Value


def Get_Sigmas_max_values(Function, Values_Result, Type_Sigmas, Sigma_x_Function, Sigma_y_Function, Tay_xy_Function,
                          z_val, W_val, U_val, V_val):
    x_array = []
    y_array = []
    z_array = [0] * Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    Max_Value = 0

    for i in range(0, Size):
        z_array[i] = [0] * Size
        x_array.append(i * step_x)
        y_array.append(i * step_y)

    for i in range(N + 1):
        if Type_Sigmas == 1:
            Sigma_x_Function = Sigma_x_Function.subs(
                [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        if Type_Sigmas == 2:
            Sigma_y_Function = Sigma_y_Function.subs(
                [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        if Type_Sigmas == 3:
            Tay_xy_Function = Tay_xy_Function.subs(
                [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])

    for j in range(0, Size):
        for i in range(0, Size):
            if Type_Sigmas == 1:
                z_array[i][j] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if Type_Sigmas == 2:
                z_array[i][j] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if Type_Sigmas == 3:
                z_array[i][j] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if abs(z_array[i][j]) > abs(Max_Value):
                Max_Value = z_array[i][j]

    return Max_Value


def Get_Sigmas_max_values_Orto(Function, Values_Result, Type_Sigmas, U_function, V_function, W_Function, z_val, W_val,
                               U_val, V_val):
    x_array = []
    y_array = []
    z_array = [0] * Size
    Max_Value = 0

    for i in range(0, Size):
        z_array[i] = [0] * Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0, Size):
        x_array.append(i * step_x)
        y_array.append(i * step_y)

    Sigma_x_Function = 5 * Xx
    Sigma_y_Function = 5 * Xx
    Tay_xy_Function = 5 * Xx
    Max_Sigmas_values = [0] * 3

    if Change == 1:
        Sigma_x_Function = (
            Get_Sigma_x_Orto(U_function, V_function, W_Function, E1_T300, nu_12_T300, nu_21_T300, z_val))
        Sigma_y_Function = (
            Get_Sigma_y_Orto(U_function, V_function, W_Function, E2_T300, nu_12_T300, nu_21_T300, z_val))
        Tay_xy_Function = (Get_Sigma_tay_Orto(U_function, V_function, W_Function, G_12_T300, z_val))
    if Change == 2:
        Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
        Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
        Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_OG, nu_OG, z_val))
    if Change == 3:
        Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))

    for i in range(N + 1):
        Sigma_x_Function = Sigma_x_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        Sigma_y_Function = Sigma_y_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])
        Tay_xy_Function = Tay_xy_Function.subs(
            [('w' + str(i), W_val[i - 1]), ('u' + str(i), U_val[i - 1]), ('v' + str(i), V_val[i - 1])])

    # print(Sigma_x_Function)
    # print(Sigma_y_Function)
    # print(Tay_xy_Function)
    Buf_1 = (1 / F1_max_T300) + (1 / F1_min_T300)
    Buf_2 = (1 / F2_max_T300) + (1 / F2_min_T300)
    Buf_3 = 1 / (F1_max_T300 * F1_min_T300)
    Buf_4 = 1 / (F2_max_T300 * F2_min_T300)
    Buf_6 = -1 / (F12_max_T300 * F12_max_T300)

    First_val = (1 / 2) * (Buf_4 - Buf_3)
    Second_Val = (1 / 2) * (Buf_3 - Buf_4)
    Third_Val = (1 / 2) * (Buf_3 + Buf_4)

    Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx, A_lenght_x / 2), (Yy, B_lenght_y / 2)])
    Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, A_lenght_x / 2), (Yy, B_lenght_y / 2)])
    Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, A_lenght_x / 2), (Yy, B_lenght_y / 2)])
    Max_Value = First_val * (Max_Sigmas_values[0] ** 2) + Second_Val * (Max_Sigmas_values[1] ** 2) + Buf_1 * \
                Max_Sigmas_values[0] + Buf_2 * Max_Sigmas_values[1] - Third_Val * (
                            (Max_Sigmas_values[0] - Max_Sigmas_values[1]) ** 2) - Buf_6 * (Max_Sigmas_values[2] ** 2)

    return Max_Value


def Get_v_coefs():
    w_coefs = []
    for i in range(1, N + 1):
        w_coefs.append(Symbol('v' + str(i)))
    return w_coefs


def Get_u_coefs():
    w_coefs = []
    for i in range(1, N + 1):
        w_coefs.append(Symbol('u' + str(i)))
    return w_coefs


def Get_w_coefs():
    w_coefs = []
    for i in range(1, N + 1):
        w_coefs.append(Symbol('w' + str(i)))
    return w_coefs


def Get_PsiX_coefs():
    w_coefs = []
    for i in range(1, N + 1):
        w_coefs.append(Symbol('PsiX' + str(i)))
    return w_coefs


def Get_PsiY_coefs():
    w_coefs = []
    for i in range(1, N + 1):
        w_coefs.append(Symbol('PsiY' + str(i)))
    return w_coefs


def Get_w_sin_x(i):
    return sin((2 * i - 1) * Xx * m.pi / A_lenght_x)


def Get_w_sin_y(j):
    return sin((2 * j - 1) * Yy * m.pi / B_lenght_y)


def Get_w_sin_x_2(i):
    return sin((2 * i) * Xx * m.pi / A_lenght_x)


def Get_w_sin_y_2(j):
    return sin((2 * j) * Yy * m.pi / B_lenght_y)


def Get_U_function_vals(u_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += u_vals[(j - 1) * N_x + i - 1] * Get_w_sin_x_2(i) * Get_w_sin_y(j)
    return Result


def Get_V_function_vals(v_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += v_vals[(j - 1) * N_x + i - 1] * Get_w_sin_x(i) * Get_w_sin_y_2(j)
    return Result


def Get_W_function_vals(w_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += w_vals[(j - 1) * N_x + i - 1] * Get_w_sin_x(i) * Get_w_sin_y(j)
    return Result


def Get_PsiX_function_vals(PsiX_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += PsiX_vals[(j - 1) * N_x + i - 1] * Get_w_sin_x(i) * Get_w_sin_y_2(j)
    return Result


def Get_PsiY_function_vals(PsiY_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += PsiY_vals[(j - 1) * N_x + i - 1] * Get_w_sin_x(i) * Get_w_sin_y(j)
    return Result


def ksi_1(PsiX_function, PsiY_function):
    if (A_Numeric):
        Result = (1 / A) * PsiX_function.diff(Xx)
    else:
        Result = (1 / A) * PsiX_function.diff(Xx) + (1 / (A * B)) * A.diff(Yy) * PsiY_function

    return Result


def ksi_2(PsiX_function, PsiY_function):
    if (A_Numeric):
        Result = (1 / B) * PsiY_function.diff(Yy)
    else:
        Result = (1 / B) * PsiY_function.diff(Yy) + (1 / (A * B)) * B.diff(Xx) * PsiX_function

    return Result


def ksi_12(PsiX_function, PsiY_function):
    Result_1 = (1 / B) * PsiX_function.diff(Yy)
    Result_2 = (1 / A) * PsiY_function.diff(Xx)
    if (A_Numeric):
        Result_3 = 0
    else:
        Result_3 = (1 / (A * B)) * (A.diff(Yy) * PsiX_function + B.diff(Xx) * PsiY_function)

    Result = (Result_1 + Result_2 - Result_3) / 2
    return Result


def Tetta_1(W_function, U_function):
    return -((1 / A) * (W_function.diff(Xx)) + K_x * U_function)


def Tetta_2(W_function, V_function):
    return -((1 / B) * (W_function.diff(Yy)) + K_y * V_function)


def Q_x(PsiX_function, PsiY_function, G_13, W_function, U_function):
    Result = k * G_13 * h * (PsiX_function - Tetta_1(W_function, U_function))
    return Result


def Q_y(PsiX_function, PsiY_function, G_23, W_function, V_function):
    Result = k * G_23 * h * (PsiY_function - Tetta_2(W_function, V_function))
    return Result


def e_x(U_function, V_function, W_Function):
    if (A_Numeric):
        Result = (1 / A) * (U_function.diff(Xx)) - K_x * W_Function + (1 / 2) * (Tetta_1(W_Function, U_function) ** 2)
    else:
        Result = (1 / A) * U_function.diff(Xx) - K_x * W_Function + (1 / 2) * (Tetta_1(W_Function, U_function) ** 2) - (
                    1 / (A * B)) * V_function * A.diff(Yy)

    return Result


def e_y(U_function, V_function, W_Function):
    if (A_Numeric):
        Result = (1 / B) * V_function.diff(Yy) - K_y * W_Function + (1 / 2) * (Tetta_2(W_Function, U_function) ** 2)
    else:
        Result = (1 / B) * V_function.diff(Yy) - K_y * W_Function + (1 / 2) * (Tetta_2(W_Function, U_function) ** 2) - (
                1 / (A * B)) * U_function * B.diff(Xx)

    return Result


def y_xy(U_function, V_function, W_Function):
    if (A_Numeric):
        Result = (1 / A) * V_function.diff(Xx) + (1 / B) * U_function.diff(Yy) + Tetta_1(W_Function,U_function) * Tetta_2(W_Function, U_function)
    else:
        Result = (1 / A) * V_function.diff(Xx) + (1 / B) * U_function.diff(Yy)
        - (1 / (A * B)) * V_function * B.diff(Xx) - (1 / (A * B)) * U_function * A.diff(Yy)
        + Tetta_1(W_Function, U_function) * Tetta_2(W_Function, U_function)

    return Result


def e_xz(z_value, U_function, V_function, W_Function, PsiX_function, PsiY_function):
    Result = e_x(U_function, V_function, W_Function) + z_value * ksi_1(PsiX_function, PsiY_function)
    return Result


def e_yz(z_value, U_function, V_function, W_Function, PsiX_function, PsiY_function):
    Result = e_y(U_function, V_function, W_Function) + z_value * ksi_2(PsiX_function, PsiY_function)
    return Result


def y_xz(z_value, U_function, V_function, W_Function, PsiX_function, PsiY_function):
    Result = y_xy(U_function, V_function, W_Function) + 2 * z_value * ksi_12(PsiX_function, PsiY_function)
    return Result


# Ortotrop functions
def N_x_Orto(z, U_function, V_function, W_Function, E_1, nu_12, nu_21):
    Result = (E_1 / (1 - nu_12 * nu_21)) * h * (
                e_x(U_function, V_function, W_Function) + nu_21 * e_y(U_function, V_function, W_Function))
    return Result


def N_y_Orto(z, U_function, V_function, W_Function, E_2, nu_12, nu_21):
    Result = (E_2 / (1 - nu_12 * nu_21)) * h * (
                e_y(U_function, V_function, W_Function) + nu_12 * e_x(U_function, V_function, W_Function))
    return Result


def N_xy_Orto(z, U_function, V_function, W_Function, G_12):
    Result = G_12 * h * y_xy(U_function, V_function, W_Function)
    return Result


def M_x_Orto(W_Function, E_1, nu_12, nu_21):
    Result = (E_1 / (1 - nu_12 * nu_21)) * ((h ** 3) / 12) * (
                ksi_1(PsiX_function, PsiY_function) + nu_21 * ksi_2(PsiX_function, PsiY_function))
    return Result


def M_y_Orto(W_Function, E_2, nu_12, nu_21):
    Result = (E_2 / (1 - nu_12 * nu_21)) * ((h ** 3) / 12) * (
                ksi_2(PsiX_function, PsiY_function) + nu_12 * ksi_1(PsiX_function, PsiY_function))
    return Result


def M_xy_Orto(W_Function, G_12):
    Result = 2 * G_12 * ((h ** 3) / 12) * ksi_12(PsiX_function, PsiY_function)
    return Result


def Get_Sigma_x_Orto(U_function, V_function, W_Function, E_1, nu_12, nu_21, z_val, PsiX_function, PsiY_function):
    Result = (E_1 / (1 - nu_12 * nu_21)) * (
                e_x(U_function, V_function, W_Function)
                + nu_21 * e_y(U_function, V_function, W_Function)
                + z_val * ksi_1(PsiX_function, PsiY_function) + nu_21 * ksi_2(PsiX_function, PsiY_function))
    return Result


def Get_Sigma_y_Orto(U_function, V_function, W_Function, E_2, nu_12, nu_21, z_val, PsiX_function, PsiY_function):
    Result = (E_2 / (1 - nu_12 * nu_21)) * (
                e_y(U_function, V_function, W_Function)
                + nu_12 * e_x(U_function, V_function, W_Function)
                + z_val * ksi_2(PsiX_function, PsiY_function) + nu_12 * ksi_1(PsiX_function, PsiY_function))
    return Result


def Get_Sigma_tay_Orto(U_function, V_function, W_Function, G_12, z_val, PsiX_function, PsiY_function):
    Result = G_12 * (y_xy(U_function, V_function, W_Function)
                     + 2 * z_val * ksi_12(PsiX_function, PsiY_function))

    return Result


# Izotrop functions
def N_x_Izo(U_function, V_function, W_Function, E, nu):
    result = (E / (1 - nu * nu)) * h * (e_x(U_function, V_function, W_Function) + nu * e_y(U_function, V_function, W_Function))
    return result


def N_y_Izo(U_function, V_function, W_Function, E, nu):
    Result = (E / (1 - nu * nu)) * h * (e_y(U_function, V_function, W_Function) + nu * e_x(U_function, V_function, W_Function))
    return Result


def N_xy_Izo(U_function, V_function, W_Function, E, nu):
    Result = (E * h * y_xy(U_function, V_function, W_Function)) / (2 * (1 + nu))
    return Result


def M_x_Izo(E, nu,PsiX_function, PsiY_function):
    Result = (E / (12 * (1 - nu * nu))) * (h ** 3) * (ksi_1(PsiX_function, PsiY_function) + nu * ksi_2(PsiX_function, PsiY_function))
    return Result


def M_y_Izo(E, nu,PsiX_function, PsiY_function):
    Result = (E / (12 * (1 - nu * nu))) * (h ** 3) * (ksi_2(PsiX_function, PsiY_function) + nu * ksi_1(PsiX_function, PsiY_function))
    return Result


def M_xy_Izo(E, nu, PsiX_function, PsiY_function):
    Result = E * ((h ** 3) / 12) * ksi_12(PsiX_function, PsiY_function) / (1 + nu)
    return Result


def Get_Sigma_x_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E, nu, z_val):
    Result = (E / (1 - nu * nu)) * (e_x(U_function, V_function, W_Function)
                                    + nu * e_y(U_function, V_function, W_Function)) \
                                    + z_val * (ksi_1(PsiX_function, PsiY_function) + nu * ksi_2(PsiX_function, PsiY_function))

    return Result


def Get_Sigma_y_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E, nu, z_val):
    Result = (E / (1 - nu * nu)) * (e_y(U_function, V_function, W_Function)
                                    + nu * e_x(U_function, V_function, W_Function)) \
                                    + z_val * (ksi_2(PsiX_function, PsiY_function) + nu * ksi_1(PsiX_function, PsiY_function))
    return Result


def Get_Sigma_tay_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E, nu, z_val):
    Result = (E / (2 + 2 * nu)) * (y_xy(U_function, V_function, W_Function)
                                    + 2 * z_val * ksi_12(PsiX_function, PsiY_function))
    return Result


def q_function(q_0, q_sv):
    Result = 0
    A1 = 0
    A_1 = [0] * 3
    A_2 = [0] * 3

    A_1[0] = 1 - ((A_lenght_x + A1) / (A_lenght_x - A1)) ** 2
    A_1[1] = 4 * (A_lenght_x + A1) / ((A_lenght_x - A1) ** 2)
    A_1[2] = -4 / ((A_lenght_x - A1) ** 2)

    A_2[0] = 1 - ((B_lenght_y + A1) / (B_lenght_y - A1)) ** 2
    A_2[1] = 4 * (B_lenght_y + A1) / ((B_lenght_y - A1) ** 2)
    A_2[2] = -4 / ((B_lenght_y - A1) ** 2)

    Result = q_0 * (A_1[0] + A_1[1] * Xx + A_1[2] * (Xx ** 2)) * (A_2[0] + A_2[1] * Yy + A_2[2] * (Yy ** 2))

    return Result


def Get_Answer(Es_Get, w_coef, u_coef, v_coef, type):

    Es = Es_Get.copy()

    if type == 2:
        Es[0] = integrate(Es[0], (Xx, Start_integral, A_lenght_x))
        Es[1] = integrate(Es[1], (Xx, Start_integral, A_lenght_x))
        Es[2] = integrate(Es[2], (Xx, Start_integral, A_lenght_x))
        Es[3] = integrate(Es[3], (Xx, Start_integral, A_lenght_x))
        Es[4] = integrate(Es[4], (Xx, Start_integral, A_lenght_x))
        Es[5] = integrate(Es[5], (Xx, Start_integral, A_lenght_x))

        Es[0] = (1 / 2) * integrate(Es[0], (Yy, Start_integral, B_lenght_y))
        Es[1] = (1 / 2) * integrate(Es[1], (Yy, Start_integral, B_lenght_y))
        Es[2] = (1 / 2) * integrate(Es[2], (Yy, Start_integral, B_lenght_y))
        Es[3] = (1 / 2) * integrate(Es[3], (Yy, Start_integral, B_lenght_y))
        Es[4] = (1 / 2) * integrate(Es[4], (Yy, Start_integral, B_lenght_y))
        Es[5] = (1 / 2) * integrate(Es[5], (Yy, Start_integral, B_lenght_y))

    Result = []
    Result_buf = []
    Buf = []
    Zeroes = [0] * N * 3
    Buf_Symbols = []

    for index in range(1, 4):
        Result.append(Result_buf)

    for i in range(1, N + 1):
        Buf.append((Es[0].diff(w_coef[i - 1]) + Es[1].diff(w_coef[i - 1]) + Es[2].diff(w_coef[i - 1]) + Es[3].diff(
            w_coef[i - 1]) + Es[4].diff(w_coef[i - 1]) + Es[5].diff(w_coef[i - 1])))
        Buf_Symbols.append(w_coef[i - 1])

    for i in range(1, N + 1):
        Buf.append((Es[0].diff(u_coef[i - 1]) + Es[1].diff(u_coef[i - 1]) + Es[2].diff(u_coef[i - 1]) + Es[3].diff(
            w_coef[i - 1]) + Es[4].diff(u_coef[i - 1]) + Es[5].diff(u_coef[i - 1])))
        Buf_Symbols.append(u_coef[i - 1])

    for i in range(1, N + 1):
        Buf.append((Es[0].diff(v_coef[i - 1]) + Es[1].diff(v_coef[i - 1]) + Es[2].diff(v_coef[i - 1]) + Es[3].diff(
            v_coef[i - 1]) + Es[4].diff(v_coef[i - 1]) + Es[5].diff(v_coef[i - 1])))
        Buf_Symbols.append(v_coef[i - 1])

    for i in range(len(Buf)):
        Buf[i] = nsimplify(Buf[i], tolerance=1e-20).evalf(15)

    for solution in linsolve(Buf, Buf_Symbols):
        Result = solution

    return Result


def Get_Jacobian(Function_E, Result_w):
    Jacobian = [0] * N * 5
    Def_Function = [0] * N * 5
    Symbol_Function = [0] * N * 5

    for i in range(0, N_x):
        Symbol_Function[i] = 'w' + str(i + 1)
        Symbol_Function[i + N_x] = 'u' + str(i + 1)
        Symbol_Function[i + N_x * 2] = 'v' + str(i + 1)
        Symbol_Function[i + N_x * 3] = 'PsiX' + str(i + 1)
        Symbol_Function[i + N_x * 4] = 'PsiY' + str(i + 1)

    for i in range(0, N * 5):
        Jacobian[i] = [0] * N * 5
        Def_Function[i] = Function_E.diff(Symbol_Function[i])
        # print("D[",i,"] =",Def_Function[i])

    for row in range(0, N * 5):
        for column in range(0, N * 5):
            Jacobian[row][column] = Def_Function[row].diff(Symbol_Function[column])
            # print("J[", row, "][",column,"] =", Jacobian[row][column])
            for W_coefs in range(0, N * 5):
                Jacobian[row][column] = Jacobian[row][column].subs((Symbol_Function[W_coefs]), (Result_w[W_coefs]))

    return Jacobian


def Get_New_iterarion(Function_E, Jackobi_inv, Last_step_X, a):
    Def_Function = [0] * N * 5

    Symbol_Function = [0] * N * 5

    for i in range(0, N_x):
        Symbol_Function[i] = 'w' + str(i + 1)
        Symbol_Function[i + N_x] = 'u' + str(i + 1)
        Symbol_Function[i + N_x * 2] = 'v' + str(i + 1)
        Symbol_Function[i + N_x * 3] = 'PsiX' + str(i + 1)
        Symbol_Function[i + N_x * 4] = 'PsiY' + str(i + 1)

    for i in range(0, N_x * 5):
        Def_Function[i] = Function_E.diff(Symbol_Function[i])
        for W_coefs in range(0, N_x * 5):
            Def_Function[i] = Def_Function[i].subs(Symbol_Function[W_coefs], (Last_step_X[W_coefs]))

    Def_Function = sym.Matrix(Def_Function)
    W = sym.Matrix(Last_step_X)
    # print("W =" ,W)
    # print("a =",a)
    # print("Ja =",Jackobi_inv)
    # print("DF =",Def_Function)
    W = W - a * (Jackobi_inv * Def_Function)
    return W


def Ne_Lin_Function_Loop(w_coefs,Es_Get,Q_Function,W_Function_get,U_function, V_function, W_Function, PsiX_function, PsiY_function):

    #symbols
    Ee = Symbol('E')
    Hh = Symbol('h')
    Qq = Symbol('q')

    print("Start loop")
    # Q_now = 0.00001
    Q_now = 0
    q_T = 3.14*2
    Q_step = q_T / 8
    Miz = 0

    q_for_graph = []
    w_for_graph = []
    w_2_for_graph = []
    W_Result = [0] * (Size + 2)
    W_Result[0] = [1] * N*5

    Es = Es_Get.copy()
    print("Start Intagrate")
    print(Es)
    #Es[0] = spi.nquad(Es[0],[Start_integral, A_lenght_x],[Start_integral,B_lenght_y])
    Es[0] = integrate(Es[0], (Xx, Start_integral, A_lenght_x))
    print("Start Intagrate 0")
    Es[1] = integrate(Es[1], (Xx, Start_integral, A_lenght_x))
    print("Start Intagrate 1")
    Es[2] = integrate(Es[2], (Xx, Start_integral, A_lenght_x))
    print("Start Intagrate 2")
    Es[3] = integrate(Es[3], (Xx, Start_integral, A_lenght_x))
    print("Start Intagrate 3")
    Es[4] = integrate(Es[4], (Xx, Start_integral, A_lenght_x))
    print("Start Intagrate 4")
    Es[5] = integrate(Es[5], (Xx, Start_integral, A_lenght_x))
    print("Start Intagrate 5")
    Es[6] = integrate(Es[6], (Xx, Start_integral, A_lenght_x))
    print("Start Intagrate 6")
    #Es[7] = integrate(Es[7], (Xx, Start_integral, A_lenght_x))

    Es[0] = (1/2) * integrate(Es[0], (Yy, Start_integral, B_lenght_y))
    print("Start Intagrate 7")
    Es[1] = (1/2) * integrate(Es[1], (Yy, Start_integral, B_lenght_y))
    print("Start Intagrate 8")
    Es[2] = (1/2) * integrate(Es[2], (Yy, Start_integral, B_lenght_y))
    print("Start Intagrate 9")
    Es[3] = (1/2) * integrate(Es[3], (Yy, Start_integral, B_lenght_y))
    print("Start Intagrate 10")
    Es[4] = (1/2) * integrate(Es[4], (Yy, Start_integral, B_lenght_y))
    print("Start Intagrate 11")
    Es[5] = (1/2) * integrate(Es[5], (Yy, Start_integral, B_lenght_y))
    print("Start Intagrate 12")
    Es[6] = (1/2) * integrate(Es[6], (Yy, Start_integral, B_lenght_y))
    print("Start Intagrate 13")
    #Es[7] = (1 / 2) * integrate(Es[7], (Yy, Start_integral, B_lenght_y))

    print("End Intagrate")
    Buf_Function = Es[0] + Es[1] + Es[2] + Es[3] + Es[4] + Es[5] + Es[6]
    Buf_Function =  Buf_Function.subs([(Ee, E), (Hh, h), (pi, m.pi)])
    print("End subs")
    j = 1
    # Цикл по разным q
    while (Miz < 1):
    #while (Q_now < 2.4):
        Q_function = q_function(Q_now, 0)
        print("Q_now = ", Q_now)
        q_for_graph.append(Q_now)

        # Es[7] = (-2) * q_function(Q_now, 0) * W_Function_get
        print("w",W_Function_get)
        Es[7] = (-2) * (Q_now * W_Function_get)
        print("Es-7 = ", Es[7])
        Es[7] = integrate(Es[7], (Xx, Start_integral, A_lenght_x))
        Es[7] = (1/2) * integrate(Es[7], (Yy, Start_integral, B_lenght_y))

        New_Buf_Function = Buf_Function + Es[7]
        print("NBF", New_Buf_Function)
        W_Result[j] = Nuton_Iter(New_Buf_Function, eps, W_Result[0], w_coefs)

        W_val = [] * N
        U_val = [] * N
        V_val = [] * N
        PsiX_val = [] * N
        PsiY_val = [] * N

        for i in range(1, N + 1):
            W_val.append(W_Result[j][i - 1])
            U_val.append(W_Result[j][i - 1 + N])
            V_val.append(W_Result[j][i - 1 + 2 * N])
            PsiX_val.append(W_Result[j][i - 1 + 3 * N])
            PsiY_val.append(W_Result[j][i - 1 + 4 * N])

        Miz = Draw_3d_Sigmas_main('Sigma_i', W_val, 1, U_function, V_function, W_Function, PsiX_function, PsiY_function,
                                  z,W_val, U_val, V_val, PsiX_val, PsiY_val)
        print("W val",W_val)
        w_for_graph.append(Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3))
        print("W = 0.5",w_for_graph[len(w_for_graph)-1])

        Q_now += Q_step
        j += 1
        print(Miz)

    Ww = [0] * 2
    Qq = [0] * 2
    Ww.append(w_for_graph[0])
    Ww.append(w_for_graph[len(w_for_graph)-1])
    Qq.append(q_for_graph[0])
    Qq.append(q_for_graph[len(q_for_graph) - 1])
    W_graph = []
    ##for i in range(1,Size+2):
    ##W_graph.append(Get_W(L/2,W_Result[i]))
    print("Plot =")
    plt.plot(w_for_graph, q_for_graph)
    plt.plot(Ww, Qq)
    plt.show()
    # axs.plot(w_for_graph,q_for_graph)
    return W_Result


def Ne_Lin_Function(w_coefs, Es_Get, Q__T):
    Ee = Symbol('E')
    Hh = Symbol('h')
    Qq = Symbol('q')
    Ll = Symbol('l')

    W_Result = [0] * N * 5

    Es = Es_Get.copy()
    Es[0] = integrate(Es[0], (Xx, Start_integral, A_lenght_x))
    Es[1] = integrate(Es[1], (Xx, Start_integral, A_lenght_x))
    Es[2] = integrate(Es[2], (Xx, Start_integral, A_lenght_x))
    Es[3] = integrate(Es[3], (Xx, Start_integral, A_lenght_x))
    Es[4] = integrate(Es[4], (Xx, Start_integral, A_lenght_x))
    Es[5] = integrate(Es[5], (Xx, Start_integral, A_lenght_x))
    Es[6] = integrate(Es[6], (Xx, Start_integral, A_lenght_x))
    Es[7] = integrate(Es[7], (Xx, Start_integral, A_lenght_x))

    Es[0] = (1 / 2) * integrate(Es[0], (Yy, Start_integral, B_lenght_y))
    Es[1] = (1 / 2) * integrate(Es[1], (Yy, Start_integral, B_lenght_y))
    Es[2] = (1 / 2) * integrate(Es[2], (Yy, Start_integral, B_lenght_y))
    Es[3] = (1 / 2) * integrate(Es[3], (Yy, Start_integral, B_lenght_y))
    Es[4] = (1 / 2) * integrate(Es[4], (Yy, Start_integral, B_lenght_y))
    Es[5] = (1 / 2) * integrate(Es[5], (Yy, Start_integral, B_lenght_y))
    Es[6] = (1 / 2) * integrate(Es[6], (Yy, Start_integral, B_lenght_y))
    Es[7] = (1 / 2) * integrate(Es[7], (Yy, Start_integral, B_lenght_y))

    Buf_Function = Es[0] + Es[1] + Es[2] + Es[3] + Es[4] + Es[5]
    Buf_Function = Buf_Function.subs([(Ee, E), (Hh, h), (pi, m.pi), (Qq, Q__T)])
    ##Main_Function = Buf_Function.subs([(Ee, E), (Hh, h), (Qq, Q__T), (Ll, 15), (pi, m.pi)])

    print(Buf_Function)

    W_Result = Nuton_Iter(Buf_Function, eps, W_Result, w_coefs)
    print("Result =", W_Result)

    return W_Result


# Nut iter
def Nuton_Iter(Function_E, eps, w0, w_coefs):
    All_Results = []
    Res_now = [0] * N * 5
    Res_Last_now = Res_now
    Res_new = []

    Now_eps = 1
    Check_Loop = 1
    Loop = 0

    Count_Iterions = 0
    New_Eps = [0] * N

    while (Now_eps > eps):
        Count_Iterions += 1
        Max_eps = 0

        Jacobi = Get_Jacobian(Function_E, Res_now)
        Jackobi_Matrix = sym.Matrix(Jacobi)
        Jacobi_Invariant = Jackobi_Matrix.inv()

        Res_new = Get_New_iterarion(Function_E, Jacobi_Invariant, Res_now, Check_Loop)
        # print("RN = ",Res_new)
        np.array(Res_new).astype(np.float64)

        for i in range(0, N):
            New_Eps[i] = abs(Res_new[i] - Res_now[i])
            if (New_Eps[i] > Max_eps):
                Max_eps = New_Eps[i]

        Res_Last_now = Res_now
        Now_eps = Max_eps
        Res_now = Res_new
        All_Results.append(Res_now)

        if Count_Iterions > 10:
            Check_Loop = Check_Loop / 10
            Res_now = Res_Last_now
            if All_Results[Count_Iterions - 3][0] < All_Results[Count_Iterions - 2][0] < \
                    All_Results[Count_Iterions - 1][0]:
                if Check_Loop != 1:
                    Check_Loop *= 10
                Res_now = Res_new
            if All_Results[Count_Iterions - 3][0] > All_Results[Count_Iterions - 2][0] > \
                    All_Results[Count_Iterions - 1][0]:
                if Check_Loop != 1:
                    Check_Loop *= 10
                Res_now = Res_new

        if Count_Iterions > 50:
            Now_eps = 0
        print("New_iteration = ", Res_now)

    return Res_now


# Main code

# Collect W-V func
w_vals = Get_w_coefs()
u_vals = Get_u_coefs()
v_vals = Get_v_coefs()
PsiX_vals = Get_PsiX_coefs()
PsiY_vals = Get_PsiY_coefs()

print("Vals")
print(w_vals)
print(u_vals)
print(v_vals)
print(PsiX_vals)
print(PsiY_vals)

W_Function = Get_W_function_vals(w_vals)
U_function = Get_U_function_vals(u_vals)
V_function = Get_V_function_vals(v_vals)
PsiX_function = Get_PsiX_function_vals(PsiX_vals)
PsiY_function = Get_PsiY_function_vals(PsiY_vals)

print("Func")
print(W_Function)
print(U_function)
print(V_function)
print(PsiX_function)
print(PsiY_function)

# Заготовочка array's
Max_W_values = [0] * 3
Max_Sigmas_values = [0] * 3
Es_main = [0] * 8
W_middle_values = []
W_middle_2_values = []
Q_values = []

# Заготовочка for results
Sigma_x_Function = 5 * Xx
Sigma_y_Function = 5 * Xx
Tay_xy_Function = 5 * Xx
Q_function = 5 * Xx

# Sigma max
QQ = 0
Check = 1

# Sigmas
U_function_buf = U_function.copy()
V_function_buf = V_function.copy()
W_function_buf = W_Function.copy()
PsiX_function_buf = PsiX_function.copy()
PsiY_function_buf = PsiY_function.copy()

z_val = -h / 2

if Change == 0:
    QQ = 0
if Change == 1:
    QQ = 1
if Change == 2:
    QQ = Sigma_t_OG
if Change == 3:
    QQ = Sigma_t_Still

# NOW
QQ = 0

if Change == 1:
    Sigma_x_Function = Get_Sigma_x_Orto(U_function_buf, V_function_buf, W_function_buf, E1_T300, nu_12_T300, nu_21_T300,
                                        z_val, PsiX_function_buf, PsiY_function_buf)
    Sigma_y_Function = Get_Sigma_y_Orto(U_function_buf, V_function_buf, W_function_buf, E2_T300, nu_12_T300, nu_21_T300,
                                        z_val, PsiX_function_buf, PsiY_function_buf)
    Tay_xy_Function = Get_Sigma_tay_Orto(U_function_buf, V_function_buf, W_function_buf, G_12_T300,
                                         z_val, PsiX_function_buf, PsiY_function_buf)
if Change == 2:
    Sigma_x_Function = Get_Sigma_x_Izo(U_function_buf, V_function_buf, W_function_buf, PsiX_function_buf, PsiY_function_buf, E1_OG, nu_OG,
                                       z_val)
    Sigma_y_Function = Get_Sigma_y_Izo(U_function_buf, V_function_buf, W_function_buf, PsiX_function_buf, PsiY_function_buf, E1_OG, nu_OG,
                                       z_val)
    Tay_xy_Function = Get_Sigma_tay_Izo(U_function_buf, V_function_buf, W_function_buf, PsiX_function_buf, PsiY_function_buf, E1_OG, nu_OG,
                                        z_val)
if Change == 3:
    Sigma_x_Function = Get_Sigma_x_Izo(U_function_buf, V_function_buf, W_function_buf, PsiX_function_buf, PsiY_function_buf, E1_Still, nu_Still,
                                       z_val)
    Sigma_y_Function = Get_Sigma_y_Izo(U_function_buf, V_function_buf, W_function_buf, PsiX_function_buf, PsiY_function_buf, E1_Still, nu_Still,
                                       z_val)
    Tay_xy_Function = Get_Sigma_tay_Izo(U_function_buf, V_function_buf, W_function_buf, PsiX_function_buf, PsiY_function_buf, E1_Still, nu_Still,
                                        z_val)

Count_num = 0
Q_now = Q_start

##dw = Lin_Function(w_coefs,y_1,q_1,Q__T)

# Es main
z_num = 0
if Change == 1:
    Es_main[0] = N_x_Orto(z_num, U_function, V_function, W_Function, E1_T300, nu_12_T300, nu_21_T300) * e_x(U_function, V_function, W_Function)
    Es_main[1] = N_y_Orto(z_num, U_function, V_function, W_Function, E2_T300, nu_12_T300, nu_21_T300) * e_y(U_function, V_function, W_Function)
    Es_main[2] = N_xy_Orto(z_num, U_function, V_function, W_Function, G_12_T300) * y_xy(U_function, V_function,W_Function)
    Es_main[3] = M_x_Orto(W_Function, E1_T300, nu_12_T300, nu_21_T300) * ksi_1(PsiX_function,PsiY_function) \
                 + M_y_Orto(W_Function,E2_T300,nu_12_T300,nu_21_T300) * ksi_2(PsiX_function, PsiY_function)
    Es_main[4] = 2 * M_xy_Orto(W_Function, G_12_T300) * ksi_12(PsiX_function, PsiY_function)
    Es_main[5] = Q_x(PsiX_function, PsiY_function, G_13_T300, W_Function, U_function)
    Es_main[6] = Q_y(PsiX_function, PsiY_function, G_23_T300, W_Function, U_function)
if Change == 2:
    Es_main[0] = N_x_Izo(U_function, V_function, W_Function, E1_OG, nu_OG) * e_x(U_function,V_function, W_Function)
    Es_main[1] = N_y_Izo(U_function, V_function, W_Function, E1_OG, nu_OG) * e_y(U_function,V_function, W_Function)
    Es_main[2] = N_xy_Izo(U_function,V_function,W_Function,E1_Still,nu_Still) * y_xy(U_function, V_function, W_Function)
    Es_main[3] = M_x_Izo(E1_OG, nu_OG,PsiX_function, PsiY_function) * ksi_1(PsiX_function, PsiY_function) \
                 + M_y_Izo(E1_OG, nu_OG,PsiX_function, PsiY_function) * ksi_2(PsiX_function, PsiY_function)
    Es_main[4] = 2 * M_xy_Izo(E1_OG, nu_OG,PsiX_function, PsiY_function) * ksi_12(PsiX_function, PsiY_function)
    Es_main[5] = Q_x(PsiX_function, PsiY_function, G_13_OG, W_Function, U_function) * (PsiX_function - Tetta_1(W_Function, U_function))
    Es_main[6] = Q_y(PsiX_function, PsiY_function, G_23_OG, W_Function, U_function) * (PsiY_function - Tetta_2(W_Function, V_function))
if Change == 3:
    Es_main[0] = N_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still) * e_x(U_function,V_function,W_Function)
    Es_main[1] = N_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still) * e_y(U_function,V_function,W_Function)
    Es_main[2] = N_xy_Izo(U_function, V_function, W_Function, E1_Still, nu_Still) * y_xy(U_function, V_function, W_Function)
    Es_main[3] = M_x_Izo(E1_Still, nu_Still,PsiX_function, PsiY_function) * ksi_1(PsiX_function, PsiY_function) \
                 + M_y_Izo(E1_Still,nu_Still,PsiX_function, PsiY_function) * ksi_2(PsiX_function, PsiY_function)
    Es_main[4] = 2 * M_xy_Izo(E1_Still, nu_Still, PsiX_function, PsiY_function) * ksi_12(PsiX_function, PsiY_function)
    #Es_main[4] = M_xy_Izo(E1_Still, nu_Still, PsiX_function, PsiY_function) * ksi_12(PsiX_function, PsiY_function)
    Es_main[5] = Q_x(PsiX_function, PsiY_function, G_13_Still, W_Function, U_function) * (PsiX_function - Tetta_1(W_Function, U_function))
    Es_main[6] = Q_y(PsiX_function, PsiY_function, G_23_Still, W_Function, V_function) * (PsiY_function - Tetta_2(W_Function, V_function))

while Check:
    Q_values.append(Q_now)
    Count_num += 1
    W_val = []
    U_val = []
    V_val = []
    PsiX_val = []
    PsiY_val = []
    W_values = []

    #Es_main[7] = (-2) * q_function(Q_now, 0) * W_Function
    Es_main_buf = Es_main.copy()
    Q_function = q_function(Q_now, 0)

    # W_values = Get_Answer(Es_main_buf, w_vals, u_vals, v_vals)
    print("Start Nelin")
    W_values = Ne_Lin_Function_Loop([0]*N*5, Es_main_buf, Q_function,W_Function,U_function, V_function, W_Function, PsiX_function, PsiY_function)
    Check = 0

    # print(W_values)
    for i in range(1, N + 1):
        W_val.append(W_values[i - 1])
        U_val.append(W_values[i - 1 + N])
        V_val.append(W_values[i - 1 + 2 * N])
        PsiX_val.append(W_values[i - 1 + 3 * N])
        PsiY_val.append(W_values[i - 1 + 4 * N])

    # for i in range(N + 1):
    # Q_function = Q_function.subs('w' + str(i), W_val[i - 1])

    # z = -h/2
    # Max_Sigmas_values[0] = Get_Sigmas_max_values('Sigma_x', W_val, 1, Sigma_x_Function, Sigma_y_Function, Tay_xy_Function, z, W_val, U_val, V_val)
    # Max_Sigmas_values[1] = Get_Sigmas_max_values('Sigma_Y', W_val, 2, Sigma_x_Function, Sigma_y_Function, Tay_xy_Function, z, W_val, U_val, V_val)
    # Max_Sigmas_values[2] = Get_Sigmas_max_values('Tay_xy',  W_val, 3, Sigma_x_Function, Sigma_y_Function, Tay_xy_Function, z, W_val, U_val, V_val)

    # Sigma_max_global=0

    # if Change == 1:
    # Sigma_max_global = abs(Get_Sigmas_max_values_Orto('Sigma_i', W_val, 1, U_function, V_function, W_Function, z, W_val, U_val, V_val))
    # print(Sigma_max_global)
    # else:
    # Sigma_max_global = ((Max_Sigmas_values[0]**2) + (Max_Sigmas_values[1]**2) - Max_Sigmas_values[0]*Max_Sigmas_values[1] + 3*(Max_Sigmas_values[2]**2))**(1/2)

    # W_middle_values.append(Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, W_Function, W_val,3))
    # W_middle_2_values.append(Get_W_Plane(A_lenght_x / 4, B_lenght_y / 4, W_Function, W_val,3))
    # Q_now += Q_max*25

    # if Sigma_max_global > QQ and Only_Value_in_Break_point != 1:
    # Check = 0
    # if Only_Value_in_Break_point and Count_num == 1:
    # Check = 0
    # Q_now -= Q_max

# Q_graph
plt.plot(W_middle_2_values, Q_values)
plt.plot(W_middle_values, Q_values)
plt.ylabel("q , МПа")
plt.xlabel('W(x,y) , м')
plt.legend(["W l/4", "W l/2"])
plt.show()
print((Q_now + Q_my) * (10 ** 6))
print(Count_num)

# Prints data :
print("W values = ", W_val)
print("U values = ", U_val)
print("V values = ", V_val)
print("Psi X values = ", V_val)
print("Psi Y values = ", V_val)

# Graph's :

Draw_3d_Q('q ', W_val, 0, Q_now, Q_function)
Max_W_values[0] = Draw_3d_W('W', W_val, 3)
Max_W_values[1] = Draw_3d_W('U', U_val, 1)
Max_W_values[2] = Draw_3d_W('V', V_val, 2)

print("max W val = ", Max_W_values[0])
print("max U val = ", Max_W_values[1])
print("max V val = ", Max_W_values[2])

z = -h / 2

if Change == 1:
    Max = Draw_3d_Sigmas_main_Orto('Критерий Хоффмана', W_val, 1, U_function, V_function, W_Function, z, W_val, U_val,
                                   V_val)
    print(abs(Max))
else:
    Draw_3d_Sigmas_main('Sigma_i', W_val, 1, U_function, V_function, W_Function, z, W_val, U_val, V_val)

Max_Sigmas_values[0] = Draw_3d_Sigmas('Sigma_x', W_val, 1, U_function, V_function, W_Function, z, W_val, U_val, V_val)
Max_Sigmas_values[1] = Draw_3d_Sigmas('Sigma_Y', W_val, 2, U_function, V_function, W_Function, z, W_val, U_val, V_val)
Max_Sigmas_values[2] = Draw_3d_Sigmas('Tay_xy', W_val, 3, U_function, V_function, W_Function, z, W_val, U_val, V_val)

print("max Sigma_x = ", Max_Sigmas_values[0])
print("max Sigma_Y = ", Max_Sigmas_values[1])
print("max Tay_xy = ", Max_Sigmas_values[2])

Sigma_max_global = ((Max_Sigmas_values[0] ** 2) + (Max_Sigmas_values[1] ** 2) - Max_Sigmas_values[0] *
                    Max_Sigmas_values[1] + 3 * (Max_Sigmas_values[2] ** 2)) ** (1 / 2)

print("Sigma max global = ", Sigma_max_global)