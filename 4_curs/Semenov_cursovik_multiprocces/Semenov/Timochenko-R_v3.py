from sympy import *
import sympy as sym
import numpy as np
import math as m
import matplotlib.pyplot as plt
from sympy import factor_terms

# Data for programm
A_Numeric = 1
Change = 3
# Data for algorithm
eps = 0.0000005
N_x = 1
N_y = 1
N = N_x * N_y
# graph points
Size = 10
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
# Symbols
Xx = Symbol('x')
Yy = Symbol('y')
# Static_DATA for variants:
# Ort:
E1 = 0
E2 = 0
nu_12 = 0
nu_21 = 0
G12 = 0
G13 = 0
G23 = 0
F1_min = 0
F1_max = 0
F2_min = 0
F2_max = 0
Density = 0
Sigma_t = 0

if Change == 1:
    # T300/976
    E1 = 1.4 * (10 ** 5)
    E2 = 0.97 * (10 ** 4)
    nu_12 = 0.29
    nu_21 = (E2 * nu_12) / E1
    G12 = 0.55 * (10 ** 4)
    G13 = G12
    G23 = 0.33 * (10 ** 4)
    F1_min = -1599
    F1_max = 1517
    F2_min = -253
    F2_max_T300 = 46
    F2_max = 41.4
    Density = 1500
    Sigma_t = 100
if Change == 2:
    # Izo:
    # Org_Glass_1:
    E1 = 2.1 * (10 ** 5)
    E2 = E1
    nu_12 = 0.30
    nu_21 = nu_12
    G12 = 0.012 * (10 ** 5)
    G13 = G12
    G23 = G12
    Density = 1190
    Sigma_t = 75
if Change == 3:
    # Still:
    E1 = 2.1 * (10 ** 5)
    E2 = E1
    nu_12 = 0.3
    nu_21 = nu_12
    #G_12_Still = 0.807 * (10 ** 5)
    G12 = E1/(2*(1+nu_12))
    G13 = G12
    G23 = G12
    Density = 7800
    Sigma_t = 265

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
def Draw_3d_Sigmas_main(Function, Values_Result, Type_Sigmas, U_function, V_function, W_Function,PsiX_function, PsiY_function, z_val, W_val, U_val,V_val, PsiX_val,PsiY_val,Sigma_x,Sigma_y,Sigma_tay):
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
        Sigma_x_Function = (Get_Sigma_x_Orto(U_function,V_function,W_Function,E1,nu_12, nu_21,z_val))
        Sigma_y_Function = (Get_Sigma_y_Orto(U_function,V_function,W_Function,E2,nu_12, nu_21,z_val))
        Tay_xy_Function = (Get_Sigma_tay_Orto(U_function,V_function,W_Function,G12,z_val))
    if Change == 2:
        Sigma_x_Function = Get_Sigma_x_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1, nu_12,z_val)
        Sigma_y_Function = Get_Sigma_y_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1, nu_12,z_val)
        Tay_xy_Function = Get_Sigma_tay_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1, nu_12,z_val)
    if Change == 3:
            #Sigma_x, Sigma_y, Sigma_tay
            #Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_Still, nu_Still, z_val))
            Sigma_x_Function = Sigma_x.copy()
            #Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_Still, nu_Still, z_val))
            Sigma_y_Function = Sigma_y.copy()
            #Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, PsiX_function, PsiY_function, E1_Still, nu_Still, z_val))
            Tay_xy_Function = Sigma_tay.copy()

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

                z_array[i][j] =((Max_Sigmas_values[0] ** 2) + (Max_Sigmas_values[1] ** 2) - Max_Sigmas_values[0] * Max_Sigmas_values[1] + 3 * (Max_Sigmas_values[2] ** 2)) ** (0.5)

                if abs(z_array[i][j]) > abs(Max_Value):
                    Max_Value = z_array[i][j]

    #plt.show()
    print("Max value Miz= ",Max_Value)
    Sigma_Tay_krit = 150
    if Change ==3:
        Sigma_Tay_krit = 340

    return Max_Value/Sigma_Tay_krit
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

def Get_w_cos_x(i):
    return cos((2 * i - 1) * Xx * m.pi / A_lenght_x)
def Get_w_cos_y(j):
    return cos((2 * j - 1) * Yy * m.pi / B_lenght_y)
def Get_w_cos_x_2(i):
    return cos((2 * i) * Xx * m.pi / A_lenght_x)
def Get_w_cos_y_2(j):
    return cos((2 * j) * Yy * m.pi / B_lenght_y)
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
            Result += PsiX_vals[(j - 1) * N_x + i - 1] * Get_w_cos_x(i) * Get_w_sin_y(j)
    return Result
def Get_PsiY_function_vals(PsiY_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += PsiY_vals[(j - 1) * N_x + i - 1] * Get_w_sin_x(i) * Get_w_cos_y(j)
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
    Result_buf = ((1 / A) * PsiY_function.diff(Xx) + (1 / B) * PsiX_function.diff(Yy))
    if (A_Numeric):
        Result_3 = 0
    else:
        Result_3 = (1 / (A * B)) * (A.diff(Yy) * PsiX_function + B.diff(Xx) * PsiY_function)

    Result = (Result_buf - Result_3) / 2
    return Result
def Tetta_1(W_function, U_function):
    Result = -((1/A) * W_function.diff(Xx)+K_x*U_function)
    return Result
def Tetta_2(W_function, V_function):
    Result = -((1/B) * W_function.diff(Yy)+K_y*V_function)
    return Result
def Q_x(PsiX_function,G_13,TETTA_1):
    Result = k * G_13 * h * (PsiX_function - TETTA_1)
    return Result
def Q_y(PsiY_function, G_23,TETTA_2):
    Result = k * G_23 * h * (PsiY_function - TETTA_2)
    return Result

def e_x(U_function, V_function, W_Function,TETTA_1):
    if (A_Numeric):
        Result = (1 / A) * U_function.diff(Xx) - K_x * W_Function + (1 / 2) * (TETTA_1*TETTA_1)
    else:
        Result = (1 / A) * U_function.diff(Xx) - K_x * W_Function + (1 / 2) * (TETTA_1 ** 2) - (
                    1 / (A * B)) * V_function * A.diff(Yy)
    return Result
def e_y(U_function, V_function, W_Function,TETTA_2):
    if (A_Numeric):
        Result = (1 / B) * V_function.diff(Yy) - K_y * W_Function + (1 / 2) * (TETTA_2 * TETTA_2)
    else:
        Result = (1 / B) * V_function.diff(Yy) - K_y * W_Function + (1 / 2) * (TETTA_2 ** 2) - (
                1 / (A * B)) * U_function * B.diff(Xx)

    return Result
def y_xy(U_function, V_function, W_Function,TETTA_1,TETTA_2):
    if (A_Numeric):
        Result = (1 / A) * V_function.diff(Xx) + (1 / B) * U_function.diff(Yy) + TETTA_1 * TETTA_2
    else:
        Result = (1 / A) * V_function.diff(Xx) + (1 / B) * U_function.diff(Yy)
        - (1 / (A * B)) * V_function * B.diff(Xx) - (1 / (A * B)) * U_function * A.diff(Yy)
        + TETTA_1 * TETTA_2

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
def N_x_Izo(E, nu,ex,ey):
    result = (E / (1 - nu * nu)) * h * (ex + nu * ey);
    return result
def N_y_Izo(E, nu,ex,ey):
    Result = (E / (1 - nu * nu)) * h * (ey + nu * ex);
    return Result
def N_xy_Izo(G, nu,yxy):
    Result = G * h * yxy
    return Result
def M_x_Izo(E, nu,ksi1,ksi2):
    Result = (E / (1 - nu * nu)) * ((h**3) / 12) * (ksi1 + nu * ksi2);
    return Result
def M_y_Izo(E, nu,ksi1,ksi2):
    Result = (E / (1 - nu * nu)) * ((h**3) / 12) * (ksi2 + nu * ksi1);
    return Result
def M_xy_Izo(G, nu, ksi12):
    Result = 2 * G * ((h**3) / 12) * ksi12
    return Result
def Get_Sigma_x_Izo(E, nu, z_val,ex,ey,KSI1,KSI2):
    Result = (E / (1 - nu * nu)) * (ex + nu * ey + z_val * (KSI1 + nu * KSI2))
    return Result
def Get_Sigma_y_Izo(E, nu, z_val,ex,ey,KSI1,KSI2):
    Result = (E / (1 - nu * nu)) * (ey + nu * ex + z_val * (KSI2 + nu * KSI1))
    return Result
def Get_Sigma_tay_Izo(G12, nu, z_val,yxy,KSI12):
    Result = G12 * (yxy + 2 * z_val * KSI12)
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
def Get_Jacobian(Function_E, Result_w):
    print("Jacodi start)")
    Jacobian = [0] * N * 5
    Def_Function = [0] * N * 5
    Symbol_Function = [0] * N * 5

    #print("Function E",Function_E)
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
            #print("J[", row, "][",column,"] =", Jacobian[row][column])
            for W_coefs in range(0, N * 5):
                Jacobian[row][column] = Jacobian[row][column].subs((Symbol_Function[W_coefs]), (Result_w[W_coefs]))
    print("Result J =", Jacobian)
    return Jacobian
def Get_New_iterarion(Function_E, Jackobi_inv, Last_step_X, a):
    Def_Function = [0] * N * 5
    Symbol_Function = [0] * N * 5
    #print( "A = ", a)
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
        #print("Def_Function[", i, "] = ", Def_Function[i])

    Def_Function = sym.Matrix(Def_Function)
    W = sym.Matrix(Last_step_X)
    #print("w = ", W)
    # print("W =" ,W)
    # print("a =",a)
    # print("Ja =",Jackobi_inv)
    # print("DF =",Def_Function)
    Result = W - a * (Jackobi_inv * Def_Function)
    #W = W * Jackobi_inv
    #print("Result = ", Result)
    return Result
def Ne_Lin_Function_Loop(w_coefs,Es_Get,W_Function_get,U_function, V_function, W_Function, PsiX_function, PsiY_function,Sigma_x,Sigma_y,Sigma_tay):
    #symbols
    Ee = Symbol('E')
    Hh = Symbol('h')
    print("Start loop")
    #Q_now = 0.01
    #Q_now = 0
    q_T = 3.14/10
    Q_step = q_T / 2
    Miz = 0
    Q_now = 0.3
    q_for_graph = []
    w_for_graph = []
    q_for_graph.append(0)
    w_for_graph.append(0)
    W_Result = [0] * (1000)
    W_Result[0] = [1] * N*5


    Es = Es_Get.copy()
    #print("while= ", Es[0])
    #Es[0] = factor_terms(Es[0])
    #print("next = ", Es[0])
    print("Start Intagrate")
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
    print(Buf_Function)
    Buf_Function =  Buf_Function.subs([(Ee, E), (Hh, h), (pi, m.pi)])
    print("End subs")

    j = 1
    # Цикл по разным q
    #while (Miz < 1):
    #for Q_now in range(1,3):
    check = 0
    while (Q_now < 3.5):
        if Q_now > 3 and check == 0:
            Q_step = Q_step/3
            check =1
        print("Q_now = ", Q_now)

        Es[7] = (-1) * (Q_now * W_Function_get)
        #print("Es-7 = ", Es[7])
        Es[7] = integrate(Es[7], (Xx, Start_integral, A_lenght_x))
        Es[7] = integrate(Es[7], (Yy, Start_integral, B_lenght_y))

        New_Buf_Function = Buf_Function + Es[7]
        W_Result[j] = Nuton_Iter(New_Buf_Function, eps, W_Result[j-1], w_coefs)

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

        print("W(",j,")",W_Result[j][i - 1])
        print("U(", j, ")", W_Result[j][i - 1 + N])
        print("V(", j, ")", W_Result[j][i - 1 + 2 * N])
        print("PsiX(", j, ")", W_Result[j][i - 1 + 3 * N])
        print("PsiY(", j, ")", W_Result[j][i - 1 + 4 * N])
        #print("Sigma_x = ",Sigma_x)
        #print("Sigma_y = ",Sigma_y)
        #print("Sigma_tay = ",Sigma_tay)
        Miz = Draw_3d_Sigmas_main('Sigma_i', W_val, 1, U_function, V_function, W_Function, PsiX_function, PsiY_function,
                                  z,W_val, U_val, V_val, PsiX_val, PsiY_val,Sigma_x,Sigma_y,Sigma_tay)

        New_w = Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3)
        Last_w = w_for_graph[len(w_for_graph)-1]
        Q_now += Q_step
        j += 1
        print("Miz = ", Miz)

        if New_w > Last_w:
            w_for_graph.append(Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3))
            q_for_graph.append(Q_now - Q_step)
        else:
            continue



    Ww = [0] * 2
    Qq = [0] * 2
    Ww.append(w_for_graph[0])
    Ww.append(w_for_graph[len(w_for_graph)-1])
    Qq.append(q_for_graph[0])
    Qq.append(q_for_graph[len(q_for_graph) - 1])
    print("Plot =")
    plt.plot(w_for_graph, q_for_graph)
    #plt.plot(Ww, Qq)
    plt.show()
    #axs.plot(w_for_graph,q_for_graph)
    return W_Result
# Nut iter
def Nuton_Iter(Function_E, eps, w0, w_coefs):
    All_Results = []
    Res_now = w0
    Now_eps = 1
    Check_Loop = 1
    Count_Iterions = 0
    New_Eps = [0] * N*5

    while (Now_eps > eps):
        Count_Iterions += 1
        Max_eps = 0
        #print("F_e", Function_E)
        Jacobi = Get_Jacobian(Function_E, Res_now)
        #print("J_e", Jacobi)
        Jackobi_Matrix = sym.Matrix(Jacobi)
        Jacobi_Invariant = Jackobi_Matrix.inv()

        #print("Jacobi = ", Jackobi_Matrix)
        Res_new = Get_New_iterarion(Function_E, Jacobi_Invariant, Res_now, Check_Loop)
        #print("Res_now = ",Res_new)
        np.array(Res_new).astype(np.float64)

        for i in range(0, N*5):
            New_Eps[i] = abs(Res_new[i] - Res_now[i])
            if (New_Eps[i] > Max_eps):
                Max_eps = New_Eps[i]

        #print("max = ",Max_eps)
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
        #print("New_iteration = ", Res_now)

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
QQ = Sigma_t
# NOW
QQ = 0
Count_num = 0
Q_now = 0

# Es main

z_num = 0
KSI_1 = ksi_1(PsiX_function,PsiY_function)
KSI_2 = ksi_2(PsiX_function,PsiY_function)
KSI_12 = ksi_12(PsiX_function,PsiY_function)
TETTA_1 = Tetta_1(W_Function,U_function)
TETTA_2 = Tetta_2(W_Function,U_function)
F_EX = e_x(U_function,V_function,W_Function,TETTA_1)
F_EY = e_y(U_function,V_function,W_Function,TETTA_2)
F_XY = y_xy(U_function,V_function,W_Function,TETTA_1, TETTA_2)


if Change == 1:
    Sigma_x_Function = Get_Sigma_x_Orto(U_function_buf, V_function_buf, W_function_buf, E1, nu_12, nu_21,z_val, PsiX_function_buf, PsiY_function_buf)
    Sigma_y_Function = Get_Sigma_y_Orto(U_function_buf, V_function_buf, W_function_buf, E2, nu_12, nu_21, z_val, PsiX_function_buf, PsiY_function_buf)
    Tay_xy_Function = Get_Sigma_tay_Orto(U_function_buf, V_function_buf, W_function_buf, G12, z_val, PsiX_function_buf, PsiY_function_buf)
if Change == 2 or Change == 3:
    Sigma_x_Function = Get_Sigma_x_Izo(E1, nu_12, z_val, F_EX, F_EY, KSI_1, KSI_2)
    Sigma_y_Function = Get_Sigma_y_Izo(E1, nu_12, z_val, F_EX, F_EY, KSI_1, KSI_2)
    Tay_xy_Function = Get_Sigma_tay_Izo(G12, nu_12, z_val, F_XY, KSI_12)

print("Sigma_x_Function = ", Sigma_x_Function)
print("Sigma_y_Function = ", Sigma_y_Function)
print("Tay_xy_Function = ", Tay_xy_Function)
integral_type = 1
if integral_type ==1:
    if Change == 1:
        Es_main[0] = N_x_Orto(z_num, U_function, V_function, W_Function, E1, nu_12, nu_21) * e_x(U_function, V_function, W_Function)
        Es_main[1] = N_y_Orto(z_num, U_function, V_function, W_Function, E1, nu_12, nu_21) * e_y(U_function, V_function, W_Function)
        Es_main[2] = N_xy_Orto(z_num, U_function, V_function, W_Function, G12) * y_xy(U_function, V_function,W_Function)
        Es_main[3] = M_x_Orto(W_Function, E1, nu_12, nu_21) * ksi_1(PsiX_function,PsiY_function) \
                     + M_y_Orto(W_Function,E2,nu_12,nu_21) * ksi_2(PsiX_function, PsiY_function)
        Es_main[4] = 2 * M_xy_Orto(W_Function, G12) * ksi_12(PsiX_function, PsiY_function)
        Es_main[5] = Q_x(PsiX_function, PsiY_function, G13, W_Function, U_function)
        Es_main[6] = Q_y(PsiX_function, PsiY_function, G23, W_Function, U_function)
    if Change == 2 or Change == 3:
        NX_I = N_x_Izo(E1, nu_12,F_EX,F_EY)
        NY_I = N_y_Izo(E1, nu_12,F_EX,F_EY)
        NXY_I = N_xy_Izo(G12, nu_12,F_XY)

        MX_I = M_x_Izo(E1, nu_12,KSI_1,KSI_2)
        MY_I = M_y_Izo(E1,nu_12,KSI_1,KSI_2)
        MXY_I = M_xy_Izo(G12, nu_12,KSI_12)

        Q_X = Q_x(PsiX_function,G13, TETTA_1)
        Q_y = Q_y(PsiY_function, G23, TETTA_2)

        Es_main[0] = NX_I * F_EX
        print("Es_main[0] = ",Es_main[0])
        Es_main[1] = NY_I * F_EY
        print("Es_main[1] = ", Es_main[1])
        Es_main[2] = NXY_I * F_XY
        print("Es_main[2] = ", Es_main[2])
        Es_main[3] = MX_I * KSI_1 + MY_I * KSI_2
        print("Es_main[3] = ", Es_main[3])
        Es_main[4] = 2 * MXY_I * KSI_12
        print("Es_main[4] = ", Es_main[4])
        Es_main[5] = Q_X * (PsiX_function - TETTA_1)
        print("Es_main[5] = ", Es_main[5])
        Es_main[6] = Q_y * (PsiY_function - TETTA_2)
        print("Es_main[6] = ", Es_main[6])
        print("New_main")

while Check:
    Q_values.append(Q_now)
    Count_num += 1
    W_val = []
    U_val = []
    V_val = []
    PsiX_val = []
    PsiY_val = []
    W_values = []

    Es_main_buf = Es_main.copy()
    Q_function = q_function(Q_now, 0)

    print("Start Nelin")
    W_values = Ne_Lin_Function_Loop([0]*N*5, Es_main_buf,W_Function,U_function, V_function, W_Function, PsiX_function, PsiY_function
                                    ,Sigma_x_Function,Sigma_y_Function,Tay_xy_Function)
    Check = 0

    for i in range(1, N + 1):
        W_val.append(W_values[i - 1])
        U_val.append(W_values[i - 1 + N])
        V_val.append(W_values[i - 1 + 2 * N])
        PsiX_val.append(W_values[i - 1 + 3 * N])
        PsiY_val.append(W_values[i - 1 + 4 * N])


    # Prints data :
    #print("W values = ", W_val)
    #print("U values = ", U_val)
    #print("V values = ", V_val)
    #print("Psi X values = ", V_val)
    #print("Psi Y values = ", V_val)
