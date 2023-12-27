import math
import os
from decimal import *
from sympy import *
import sympy as sym
import numpy as np
import math as m
import matplotlib.pyplot as plt
from sympy import factor_terms
import time
import multiprocessing as mp
import mpmath as vpa

#m to mm
# Data for programm
One_point = 0
Multi = 1
A_Numeric = 1
Change = 2
Type_Resolve = 1# 0 - Ritz 1 - Nuton 2 - Bubnov
# Data for algorithm
Start_w0 = 0.0035
eps = 0.05
N_x = 2
N_y = 2
N = N_x * N_y
q_find = 0.6/3
Q_Count_Steps = 30
Q_step = q_find/Q_Count_Steps
Count_Es = 12
# graph points
Size = 30
Size_visualization = 10
# Data about object main
H_coef = 1
h = 0.22 * H_coef  #mm
A_lenght_x = 200
B_lenght_y = 200
R_1_x = 5000
R_2_y = 3330
A = 1
B = 1
K_x = 1 / R_1_x
K_y = 1 / R_2_y
z = 0
# Data material
V = 200 * 200 * 0.22  # m^3
P1 = 1500
P2 = 7800
P3 = 1190
Square = 20 * 20
k = 5 / 6
# Symbols
Xx = Symbol('x')
Yy = Symbol('y')
Hh = Symbol('h')
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
    E1 = 3000000
    E2 = E1
    nu_12 = 0.35
    nu_21 = nu_12
    G12 = 1200000
    G13 = G12
    G23 = G12
    Density = 1190
    Sigma_t = 75000
if Change == 3:
    # Still:
    E1 = 210000000
    E2 = E1
    nu_12 = 0.3
    nu_21 = nu_12
    G12 = 80700000
    G13 = G12
    G23 = G12
    Density = 7800
    Sigma_t = 265000

# Settings to integral
Start_integral = 0
def intagrete_Es(fun,queue,i):
    timer = time.time()
    print("Start Intagrate = ",i)
    Result = fun
    Result = factor_terms(Result)
    Result = nsimplify(Result, tolerance=1e-15).evalf(15)
    Result = integrate(Result, (Xx, Start_integral, A_lenght_x))
    Result = factor_terms(Result)
    Result = nsimplify(Result, tolerance=1e-15).evalf(15)
    Result = (1 / 2) * integrate(Result, (Yy, Start_integral, B_lenght_y))
    Result = nsimplify(Result, tolerance=1e-15).evalf(15)
    queue.put(Result)
    spend_time = time.time() - timer
    print("Result[",i," By time = ",spend_time,"=",Result,)
    return Result
def intagrete_Es_2(data):
    fun = data[0]
    i = data[1]
    #timer = time.time()
    #print("Start Intagrate = ",i)
    Result = fun
    Result = factor_terms(Result)
    #Result = nsimplify(Result, tolerance=1e-15).evalf(15)
    Result = integrate(Result, (Xx, Start_integral, A_lenght_x))
    Result = factor_terms(Result)
    #Result = nsimplify(Result, tolerance=1e-15).evalf(15)
    Result = (1 / 2) * integrate(Result, (Yy, Start_integral, B_lenght_y))
    #Result = nsimplify(Result, tolerance=1e-15).evalf(15)
    #queue.put(Result)
    #spend_time = time.time() - timer
    #print("Result[",i," By time = ",spend_time,"=",Result,)

    Result = Result.expand()
    timer = time.time()
    Result = factor_terms(Result)

    return Result
def intagrete_Es_Simp(fun,queue,i):
    timer = time.time()
    print("Start Intagrate = ", i)
    Result = fun
    #Result = factor_terms(Result)

    Intervals_count = 8

    h1 = A_lenght_x/Intervals_count
    h2 = B_lenght_y/Intervals_count

    integrate_x = 0
    integrate_x += Result.subs(Xx,0)
    integrate_x += Result.subs(Xx,A_lenght_x)
    integrate_x = integrate_x.evalf()

    for Index in range(1, Intervals_count):
        if Index % 2:
            integrate_x += 4 * Result.subs(Xx, Index * h1)
        else:
            integrate_x += 2 * Result.subs(Xx, Index * h1)
    integrate_x *= h1/3
    Buf_function = integrate_x.evalf()

    #Buf_function =
    #Result = integrate_x
    #Result = nsimplify(integrate_x, tolerance=1e-15).evalf(15)
    Buf_function = factor_terms(Buf_function)
    integrate_y = 0
    integrate_y += Buf_function.subs(Yy, 0)
    integrate_y += Buf_function.subs(Yy, B_lenght_y)
    integrate_y = integrate_y.evalf()

    for Index in range(1,Intervals_count):
        if Index % 2:
            integrate_y += 4 * Buf_function.subs(Yy, Index * h2)
        else:
            integrate_y += 2 * Buf_function.subs(Yy, Index * h2)

    integrate_y *= h2/3

    Result = integrate_y.evalf()
    #Result = factor_terms(Result)
    #Result = nsimplify(Result, tolerance=1e-10).evalf(10)
    if N == 1:
        Result = Result.simplify()
    getcontext().prec = 6
    #Result = Decimal(Result)
    #Result = nsimplify(Result, tolerance=1e-22).evalf(22)
    queue.put(Result)
    spend_time = time.time() - timer
    #print("Result[",i," By time = ",spend_time,"=",Result,)
    print("Result[", i, " By time = ", spend_time )
def intagrete_Es_Simp_2(data):
    fun = data[0]
    i = data[1]
    timer = time.time()
    #print("Start Intagrate = ", i)
    Result = fun
    #Result = factor_terms(Result)

    Intervals_count = 5

    h1 = A_lenght_x/Intervals_count
    h2 = B_lenght_y/Intervals_count

    integrate_x = 0
    integrate_x += Result.subs(Xx,0)
    integrate_x += Result.subs(Xx,A_lenght_x)
    integrate_x = integrate_x.evalf(15)

    for Index in range(1, Intervals_count):
        if Index % 2:
            integrate_x += 4 * Result.subs(Xx, Index * h1)
        else:
            integrate_x += 2 * Result.subs(Xx, Index * h1)
    integrate_x *= h1/3
    Buf_function = integrate_x.evalf(15)

    #Buf_function =
    #Result = integrate_x
    #Result = nsimplify(integrate_x, tolerance=1e-15).evalf(15)
    #Buf_function = factor_terms(Buf_function)
    integrate_y = 0
    integrate_y += Buf_function.subs(Yy, 0)
    integrate_y += Buf_function.subs(Yy, B_lenght_y)
    integrate_y = integrate_y.evalf(15)

    for Index in range(1,Intervals_count):
        if Index % 2:
            integrate_y += 4 * Buf_function.subs(Yy, Index * h2)
        else:
            integrate_y += 2 * Buf_function.subs(Yy, Index * h2)

    integrate_y *= h2/3

    Result = integrate_y.evalf(15)
    #Result = factor_terms(Result)
    #Result = nsimplify(Result, tolerance=1e-10).evalf(10)
    #if N == 1:
        #Result = Result.simplify()
    #Result = Result.simplify()
    #mp.dps = 5
    vpa.mp.dps = 8
    #Result = vpa.mpf(Result)
    #Result = trigsimp(Result)
    #Result = Decimal(Result)
    #Result = nsimplify(Result, tolerance=1e-22).evalf(22)
    #queue.put(Result)
    Result = Result.expand()
    #timer = time.time()
    Result = factor_terms(Result)
    #Result = Result.simplify()
    #print("Simplify =  By time = ", time.time() - timer )
    #spend_time = time.time() - timer
    #print("Result[",i," By time = ",spend_time,"=",Result,)
    #print("Result[", i, "] By time = ", spend_time )
    return Result
def get_diff(fun,queue,symbol,i):
    timer = time.time()
    print("start diff =",i)
    Result = []
    Result.append(fun.diff(symbol))
    Result.append(i)
    queue.put(Result)
    spend_time = time.time() - timer
    print("Result[", i, " By time = ", spend_time)
    return 0
def get_diff_2(data):
    i = data[2]
    symbol = data[1]
    fun = data[0]
    timer = time.time()
    #print("start diff =",i)
    Result = []
    Result = fun.diff(symbol)
    #Result.append(i)
    #queue.put(Result)
    spend_time = time.time() - timer
    #print("Result[", i, " By time = ", spend_time)
    return Result
def get_diff_2_Ja(data):
    i = data[2]
    symbol = data[1]
    fun = data[0]
    timer = time.time()
    #print("start diff =",i)
    Result = []
    buf_result = 0
    for index in range(0,Count_Es):
        buf_result+=fun[index].diff(symbol)
    Result = buf_result
    #Result.append(i)
    #queue.put(Result)
    spend_time = time.time() - timer
    #print("Result[", i, " By time = ", spend_time)
    return Result

def get_diff_2arg(fun,queue,symbol,i,j):
    timer = time.time()
    print("start diff =",i,",",j)
    Result = []
    Result.append(fun.diff(symbol))
    Result.append(i)
    Result.append(j)
    queue.put(Result)
    spend_time = time.time() - timer
    print("Result[", i,",",j, " By time = ", spend_time)
    return 0
def get_diff_2arg_2(data):
    fun = data[0]
    symbol = data[1]
    i = data[2]
    j = data[3]
    #print(i, j)
    timer = time.time()
    #print("start diff =",i,",",j)
    Result = []
    Result.append(fun.diff(symbol))
    Result.append((i,j))
    #Result.append(j)
    #queue.put(Result)
    spend_time = time.time() - timer
    #print("Result[", i,",",j, " By time = ", spend_time)
    return Result
def get_diff_2arg_2_Ja(data):
    fun = data[0]
    symbol = data[1]
    i = data[2]
    j = data[3]
    #print(i, j)
    timer = time.time()
    #print("start diff =",i,",",j)
    Result = []
    buf_result = 0
    for index in range(0,Count_Es):
        buf_result+=fun[index].diff(symbol)
    Result.append(buf_result)
    Result.append((i,j))
    #Result.append(j)
    #queue.put(Result)
    spend_time = time.time() - timer
    #print("Result[", i,",",j, " By time = ", spend_time)
    return Result

def Sub_data(data):
    fun = data[0]
    Symbols_get = data[1]
    #print("Sym_data =",Symbols_get)
    Result = fun.subs(Symbols_get)
    return Result
def intagrete_Es_Simp_Bubnov(fun, queue, i):
    timer = time.time()
    print("Start Intagrate = ", i)
    Result = fun
    Result = factor_terms(Result)

    Intervals_count = 10

    h1 = A_lenght_x / Intervals_count
    h2 = B_lenght_y / Intervals_count

    integrate_x = 0
    integrate_x += Result.subs(Xx, 0)
    integrate_x += Result.subs(Xx, A_lenght_x)
    integrate_x = integrate_x.evalf()

    for Index in range(1, Intervals_count):
        if Index % 2:
            integrate_x += 4 * Result.subs(Xx, Index * h1)
        else:
            integrate_x += 2 * Result.subs(Xx, Index * h1)
    integrate_x *= h1 / 3
    Buf_function = integrate_x.evalf()

    # Buf_function =
    # Result = integrate_x
    # Result = nsimplify(integrate_x, tolerance=1e-15).evalf(15)
    Buf_function = factor_terms(Buf_function)
    integrate_y = 0
    integrate_y += Buf_function.subs(Yy, 0)
    integrate_y += Buf_function.subs(Yy, B_lenght_y)
    integrate_y = integrate_y.evalf()

    for Index in range(1, Intervals_count):
        if Index % 2:
            integrate_y += 4 * Buf_function.subs(Yy, Index * h2)
        else:
            integrate_y += 2 * Buf_function.subs(Yy, Index * h2)

    integrate_y *= h2 / 3

    Result = integrate_y.evalf()
    # Result = factor_terms(Result)
    # Result = nsimplify(Result, tolerance=1e-15).evalf(15)
    Result = Result.simplify()
    Result = nsimplify(Result, tolerance=1e-22).evalf(22)
    return_val = [0]*2
    return_val[0] = Result
    return_val[1] = i
    queue.put(return_val)
    spend_time = time.time() - timer
    # print("Result[",i," By time = ",spend_time,"=",Result,)
    print("Result[", i, " By time = ", spend_time)
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
def Draw_deformed_plot(W_func,U_func,V_func):
    x_array = [0] * (Size_visualization*Size_visualization)
    y_array = [0] * (Size_visualization*Size_visualization)

    x_array_graph = [0] * (Size_visualization*Size_visualization)
    y_array_graph = [0] * (Size_visualization*Size_visualization)
    pi = 3.14

    Min_ab = min(A_lenght_x,B_lenght_y)
    Max_ab = max(A_lenght_x,B_lenght_y)

    Min_ab_g = min(pi * R_1_x, pi * R_2_y)
    Max_ab_g = max(pi * R_1_x, pi * R_2_y)

    X_array = [0] * (Size_visualization * Size_visualization)
    Y_array = [0] * (Size_visualization * Size_visualization)
    Z_array = [] * (Size_visualization * Size_visualization)

    R_main = max(R_1_x,R_2_y) - min(R_1_x,R_2_y)
    r_main = min(R_1_x,R_2_y)

    x_step = (2*Min_ab / (2 * r_main)) / (Size_visualization**2)
    y_step = (2*Max_ab / (2 * (R_main + r_main))) / (Size_visualization**2)
    print("X_step =",x_step)
    print("Y_step = ",y_step)

    x_array[0] = -Min_ab/(2*r_main)
    y_array[0] = -Max_ab/(2*(R_main+r_main))
    x_array_graph[0] = 0
    y_array_graph[0] = 0
    x_graph_step = (2*Min_ab_g / (2 * r_main)) / (Size_visualization**2)
    y_graph_step = (2*Max_ab_g / (2 * (R_main + r_main))) / (Size_visualization ** 2)

    for index in range(1,Size_visualization*Size_visualization):
        x_array[index] = x_array[index - 1] + x_step
        y_array[index] = y_array[index - 1] + y_step
    for index in range(1,Size_visualization*Size_visualization):
        x_array_graph[index] = x_array_graph[index - 1] + x_graph_step
        y_array_graph[index] = y_array_graph[index - 1] + y_graph_step

    print(x_array)
    print(y_array)
    for index in range(0,Size_visualization**2):
        buf = [0] * (Size_visualization**2)
        Z_array.append(buf)
        for index_2 in range(0,Size_visualization**2):
            Z_array[index][index_2] = ((r_main*sin(x_array_graph[index_2]) + r_main*sin(y_array_graph[index]))/40)*h

    print("Z0 = " ,Z_array[0])
    print("Z_len", len(Z_array))
    Z = np.array(Z_array)
    x_array_graph = np.array(x_array_graph)
    y_array_graph = np.array(y_array_graph)
    X, Y = np.meshgrid(x_array, y_array)

    #print("X[0][0]",X[0][0])
    #for i in range(0,len(X)):
        #print(X[i])

    fig = plt.figure()
    ax = plt.axes(projection='3d',ylim =(-0.02,0.02),xlim = (-0.02,0.02),zlim = (-0.02,0.02))
    ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x [м]')
    ax.set_ylabel('y [м]')
    ax.set_zlabel('z [м]')
    ax.set_title("Visual")
    plt.show()

    Nu_array = [0] * 3
    Nv_array = [0] * 3
    Nw_array = [0] * 3

    Nu_array[0] = -sin(Xx)*sin(Yy)*U_func
    Nu_array[1] = -sin(Xx)*cos(Yy)*U_func
    Nu_array[2] = cos(Xx)*V_func

    Nv_array[0] = cos(Yy) * V_func
    Nv_array[1] = -sin(Yy) * V_func
    Nv_array[2] = 0 * V_func

    Nw_array[0] = cos(Xx)*sin(Yy) * W_func
    Nw_array[1] = cos(Xx)*cos(Yy) * W_func
    Nw_array[2] = sin(Xx) * W_func

    print(W_func)
    print(U_func)
    print(V_func)
    print(x_array)
    x_new = np.linspace(0,A_lenght_x,Size_visualization*Size_visualization)
    y_new = np.linspace(0,B_lenght_y,Size_visualization*Size_visualization)
    for index in range(0,Size_visualization*Size_visualization):
        for index_2 in range(0, Size_visualization * Size_visualization):
            #X[index][index_2] = X[index][index_2] + Nu_array[0].subs([(Xx, x_array[index_2]), (Yy, y_array[index])]) \
            #                          + Nv_array[0].subs([(Xx, x_array[index_2]), (Yy, y_array[index])]) + Nw_array[0].subs([(Xx, x_array[index_2]), (Yy, y_array[index])])
            #Y[index][index_2] = Y[index][index_2] + Nu_array[1].subs([(Xx, x_array[index_2]), (Yy, y_array[index])]) \
            #                          + Nv_array[2].subs([(Xx, x_array[index_2]), (Yy, y_array[index])]) + Nw_array[2].subs([(Xx, x_array[index_2]), (Yy, y_array[index])])
            #Z_array[index][index_2] = Z_array[index][index_2]/10 + Nw_array[1].subs([(Xx, x_new[index_2]), (Yy, y_new[index])]) \
            #                          + Nw_array[1].subs([(Xx, x_new[index_2]), (Yy, y_new[index])]) + Nw_array[2].subs([(Xx, x_new[index_2]), (Yy, y_new[index])])
            X[index][index_2] = X[index][index_2] - U_func.subs([(Xx, x_new[index]), (Yy, y_new[index_2])])*10
            Y[index][index_2] = Y[index][index_2] - V_func.subs( [(Xx, x_new[index]), (Yy, y_new[index_2])])*10
            Z_array[index][index_2] =Z_array[index][index_2] - W_func.subs( [(Xx, x_new[index]), (Yy, y_new[index_2])])*50
            #print("old",Z_array[index][index_2]/1000)
            #print("new",Nu_array[2].subs([(Xx, x_array[index_2]), (Yy, y_array[index])]) \
            #                          + Nv_array[2].subs([(Xx, x_array[index_2]), (Yy, y_array[index])]) + Nw_array[2].subs([(Xx, x_array[index_2]), (Yy, y_array[index])]))

    #z_array = np.array(d_z_array)
    #x_array_graph = np.array(x_array_graph)
    #y_array_graph = np.array(y_array_graph)
    #X, Y = np.meshgrid(x_array_graph, y_array_graph)
    Z = np.array(Z_array)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x [м]')
    ax.set_ylabel('y [м]')
    ax.set_zlabel('z [м]')
    ax.set_title("Visual")
    plt.show()
def Draw_3d_W(Function, Values_Result, type_f):
    x_array = []
    y_array = []
    z_array = [0] * Size
    Max_Value = 0
    Get_W = Function.copy()
    # print(Values_Result)
    print(Function)
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
            z_array[i][j] = Get_W.subs( [(Xx, x_array[i]), (Yy, y_array[j])])

    z_array = np.array(z_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x [м]')
    ax.set_ylabel('y [м]')
    ax.set_zlabel("W" + ' (x,y) [м]')
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
def get_Mizz(W_val, U_val,V_val, PsiX_val,PsiY_val,Sigma_x,Sigma_y,Sigma_tay):
    x_array = []
    y_array = []
    z_array = [0] * Size

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
        print("Change == 1")
    if Change == 2:
        print("Change == 2")
    if Change == 3:
            Sigma_x_Function = Sigma_x.copy()
            Sigma_y_Function = Sigma_y.copy()
            Tay_xy_Function = Sigma_tay.copy()

    #timer = time.time()
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
    #print("Spend time to Mizz =",time.time() - timer)
    x_Mizz = 0
    y_Mizz = 0
    Max_Value = 0
    Mx_max = 0
    My_max = 0
    mt_max = 0
    for j in range(0, Size):
        for i in range(0, Size):
                Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                if abs(Max_Sigmas_values[0])>Mx_max:
                    Mx_max = abs(Max_Sigmas_values[0])
                if abs(Max_Sigmas_values[1])>My_max:
                    My_max = abs(Max_Sigmas_values[1])
                if abs(Max_Sigmas_values[2])>mt_max:
                    mt_max = abs(Max_Sigmas_values[2])

                z_array[i][j] =((Max_Sigmas_values[0] ** 2) + (Max_Sigmas_values[1] ** 2) - Max_Sigmas_values[0] * Max_Sigmas_values[1] + 3*(Max_Sigmas_values[2] ** 2)) ** (0.5)
                #if i == 10 and j == 10:
                    #print("Mizz in center =",z_array[i][j])
                if abs(z_array[i][j]) > abs(Max_Value):
                    Max_Value = z_array[i][j]
                    x_Mizz = x_array[i]
                    y_Mizz = y_array[j]

    print("Mx  =", Mx_max)
    print("My  =", My_max)
    print("Mtay  =", mt_max)
    #plt.show()
    #print("Max value Miz= ",Max_Value)
    #Sigma_Tay_krit = 150
    #if Change ==3:
        #Sigma_t = 300

    return Max_Value,x_Mizz,y_Mizz
def get_Mizz_2(W_val, U_val,V_val, PsiX_val,PsiY_val,Sigma_x,Sigma_y,Sigma_tay):
    x_array = []
    y_array = []
    z_array = [0] * Size

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
        print("Change == 1")
    if Change == 2:
        print("Change == 2")
    if Change == 3:
            Sigma_x_Function = Sigma_x.copy()
            Sigma_y_Function = Sigma_y.copy()
            Tay_xy_Function = Sigma_tay.copy()

    #timer = time.time()
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
    #print("Spend time to Mizz =",time.time() - timer)
    x_Mizz = 0
    y_Mizz = 0
    Max_Value = 0
    Mx_max = 0
    My_max = 0
    mt_max = 0
    for j in range(0, Size):
        for i in range(0, Size):
                Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                if abs(Max_Sigmas_values[0])>Mx_max:
                    Mx_max = abs(Max_Sigmas_values[0])
                if abs(Max_Sigmas_values[1])>My_max:
                    My_max = abs(Max_Sigmas_values[1])
                #print(mt_max)
                #print(Max_Sigmas_values[2])
                if abs(Max_Sigmas_values[2])>mt_max:
                    mt_max = abs(Max_Sigmas_values[2])

                z_array[i][j] =((Max_Sigmas_values[0] - Max_Sigmas_values[1] )**2 + (Max_Sigmas_values[1] - 0) ** 2 + (Max_Sigmas_values[0]) ** 2 + 6*(Max_Sigmas_values[2] ** 2)) ** (0.5)
                z_array[i][j] = z_array[i][j] * 1/(2 **(0.5))
                #if i == 10 and j == 10:
                    #print("Mizz in center =",z_array[i][j])
                if z_array[i][j] > Max_Value:
                    Max_Value = z_array[i][j]
                    x_Mizz = x_array[i]
                    y_Mizz = y_array[j]

    print("Mx  =", Mx_max)
    print("My  =", My_max)
    print("Mtay  =", mt_max)
    #plt.show()
    #print("Max value Miz= ",Max_Value)
    #Sigma_Tay_krit = 150
    #if Change ==3:
        #Sigma_t = 300

    return Max_Value/Sigma_t,x_Mizz,y_Mizz

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
        Result = PsiX_function.diff(Xx)/A
    else:
        Result = (1 / A) * PsiX_function.diff(Xx) + (1 / (A * B)) * A.diff(Yy) * PsiY_function

    return Result
def ksi_2(PsiX_function, PsiY_function):
    if (A_Numeric):
        Result = PsiY_function.diff(Yy)/B
    else:
        Result = (1 / B) * PsiY_function.diff(Yy) + (1 / (A * B)) * B.diff(Xx) * PsiX_function

    return Result
def ksi_12(PsiX_function, PsiY_function):
    Result_buf = (PsiY_function.diff(Xx)/A + PsiX_function.diff(Yy)/B)
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
    Result = k * G_13 * Hh * (PsiX_function - TETTA_1)
    return Result
def Q_y(PsiY_function, G_23,TETTA_2):
    Result = k * G_23 * Hh * (PsiY_function - TETTA_2)
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
def e_x_z(ex,ksi_1):
    Result = ex + (-Hh/2)*ksi_1
    return Result
def e_y_z(ey,ksi_2):
    Result = ey + (-Hh/2)*ksi_2
    return Result
def y_xy_z(yxy,ksi_12):
    Result = yxy + (-Hh)*ksi_12

    return Result

def y_xy_part(U_function, V_function, W_Function,TETTA_1,TETTA_2):
    Result = [0] * 3
    Result[0] = (1 / A) * V_function.diff(Xx)
    Result[1] = (1 / B) * U_function.diff(Yy)
    Result[2] = TETTA_1 * TETTA_2

    return Result
# Ortotrop functions
def N_x_Orto(z, U_function, V_function, W_Function, E_1, nu_12, nu_21):
    Result = (E_1 / (1 - nu_12 * nu_21)) * Hh * (
                e_x(U_function, V_function, W_Function) + nu_21 * e_y(U_function, V_function, W_Function))
    return Result
def N_y_Orto(z, U_function, V_function, W_Function, E_2, nu_12, nu_21):
    Result = (E_2 / (1 - nu_12 * nu_21)) * Hh * (
                e_y(U_function, V_function, W_Function) + nu_12 * e_x(U_function, V_function, W_Function))
    return Result
def N_xy_Orto(z, U_function, V_function, W_Function, G_12):
    Result = G_12 * Hh * y_xy(U_function, V_function, W_Function)
    return Result
def M_x_Orto(W_Function, E_1, nu_12, nu_21):
    Result = (E_1 / (1 - nu_12 * nu_21)) * ((Hh ** 3) / 12) * (
                ksi_1(PsiX_function, PsiY_function) + nu_21 * ksi_2(PsiX_function, PsiY_function))
    return Result
def M_y_Orto(W_Function, E_2, nu_12, nu_21):
    Result = (E_2 / (1 - nu_12 * nu_21)) * ((Hh ** 3) / 12) * (
                ksi_2(PsiX_function, PsiY_function) + nu_12 * ksi_1(PsiX_function, PsiY_function))
    return Result
def M_xy_Orto(W_Function, G_12):
    Result = 2 * G_12 * ((Hh ** 3) / 12) * ksi_12(PsiX_function, PsiY_function)
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
def N_x_Izo(E, nu12,nu21,ex,ey):
    result = (E / (1 - nu12 * nu21)) * Hh * (ex + nu12 * ey);
    return result
def N_x_Izo_part(E, nu12,nu21,ex,ey):
    Result = [0]*2
    Result[0] = (E / (1 - nu12 * nu21)) * Hh * ex;
    Result[1] = (E / (1 - nu12 * nu21)) * Hh * nu12 * ey;
    return Result
def N_y_Izo(E, nu12,nu21,ex,ey):
    Result = (E / (1 - nu12 * nu21)) * Hh * (ey + nu21 * ex);
    return Result
def N_y_Izo_part(E, nu12,nu21,ex,ey):
    Result = [0]*2
    Result[0] = (E / (1 - nu12 * nu21)) * Hh * ey;
    Result[1] = (E / (1 - nu12 * nu21)) * Hh * nu21 * ex;
    return Result
def N_xy_Izo(G, nu,yxy):
    Result = [0]*3
    Result[0] = G * Hh * yxy[0]
    Result[1] = G * Hh * yxy[1]
    Result[2] = G * Hh * yxy[2]
    return Result
def M_x_Izo(E, nu,ksi1,ksi2):
    Result = (E / (1 - nu * nu)) * ((Hh**3) / 12) * (ksi1 + nu * ksi2);
    return Result
def M_y_Izo(E, nu,ksi1,ksi2):
    Result = (E / (1 - nu * nu)) * ((Hh**3) / 12) * (ksi2 + nu * ksi1);
    return Result
def M_xy_Izo(G, nu, ksi12):
    Result = 2 * G * ((Hh**3) / 12) * ksi12
    return Result
def Get_Sigma_x_Izo(E, nu, z_val,ex,ey,KSI1,KSI2):
    Result = (E / (1 - nu * nu)) * (ex + nu * ey + (-Hh/2) * (KSI1 + nu * KSI2))
    return Result
def Get_Sigma_y_Izo(E, nu, z_val,ex,ey,KSI1,KSI2):
    Result = (E / (1 - nu * nu)) * (ey + nu * ex + (-Hh/2) * (KSI2 + nu * KSI1))
    return Result
def Get_Sigma_tay_Izo(G12, nu, z_val,yxy,KSI12):
    Result = G12 * (yxy + 2 * (-Hh/2) * KSI12)
    return Result
def q_function(q_0, q_sv):
    Qq = Symbol('q')
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

    Result = Qq * (A_1[0] + A_1[1] * Xx + A_1[2] * (Xx ** 2)) * (A_2[0] + A_2[1] * Yy + A_2[2] * (Yy ** 2))

    return Result
#Rits
def Get_Answer_Rits(Es_Get,w_coef,u_coef,v_coef):
    Es = Es_Get.copy()

    Es[0] = integrate(Es[0], (Xx, Start_integral, A_lenght_x))
    Es[1] = integrate(Es[1], (Xx, Start_integral, A_lenght_x))
    Es[2] = integrate(Es[2], (Xx, Start_integral, A_lenght_x))
    Es[3] = integrate(Es[3], (Xx, Start_integral, A_lenght_x))
    Es[4] = integrate(Es[4], (Xx, Start_integral, A_lenght_x))
    Es[5] = integrate(Es[5], (Xx, Start_integral, A_lenght_x))

    Es[0] = (1/2) * integrate(Es[0], (Yy, Start_integral, B_lenght_y))
    Es[1] = (1/2) * integrate(Es[1], (Yy, Start_integral, B_lenght_y))
    Es[2] = (1/2) * integrate(Es[2], (Yy, Start_integral, B_lenght_y))
    Es[3] = (1/2) * integrate(Es[3], (Yy, Start_integral, B_lenght_y))
    Es[4] = (1/2) * integrate(Es[4], (Yy, Start_integral, B_lenght_y))
    Es[5] = (1/2) * integrate(Es[5], (Yy, Start_integral, B_lenght_y))

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
#Bub-Gal
def get_L1(Nx,Ny,Nxy,Qx,Tetta_1,Tetta_2):
    if (A_Numeric):
        Result = Nx.diff(Xx) + Nxy.diff(Yy) - K_x*Qx +K_x * (Nx*Tetta_1 + Nxy * Tetta_2)
    else:
        Result = 0
    return Result
def get_L2(Nx,Ny,Nxy,Qy,Tetta_1,Tetta_2):
    if (A_Numeric):
        Result = Ny.diff(Yy) + Nxy.diff(Xx) - K_y*Qy + K_y * (Ny*Tetta_2 + Nxy * Tetta_1)
    else:
        Result = 0
    return Result
def get_L3(Nx,Ny,Nxy,Tetta_1,Tetta_2,Qx,Qy,q):
    if (A_Numeric):
        Result = K_x*Nx + K_y*Ny - (Nx*Tetta_1 + Nxy*Tetta_2).diff(Xx) \
                 - (Ny*Tetta_2 + Nxy*Tetta_1).diff(Yy) + Qx.diff(Xx) + Qy.diff(Yy) + q
    else:
        Result = 0
    return Result
def get_L4(Mx,Mxy,Qx):
    if (A_Numeric):
        Result = Mx.diff(Xx) + Mxy.diff(Yy) - Qx

    else:
        Result = 0
    return Result
def get_L5(My,Mxy,Qy):
    if (A_Numeric):
        Result = My.diff(Yy) + Mxy.diff(Xx) - Qy
    else:
        Result = 0
    return Result
def Get_Es_Bub(Nx,Ny,Nxy,Mx,My,Mxy,Qx,Qy,Tetta_1,Tetta_2,Q_find):
    Es = [0]*N*5

    L1 = get_L1(Nx,Ny,Nxy,Qx,Tetta_1,Tetta_2)
    L2 = get_L2(Nx,Ny,Nxy,Qy,Tetta_1,Tetta_2)
    L3 = get_L3(Nx, Ny, Nxy, Tetta_1, Tetta_2, Qx, Qy, 0)
    L4 = get_L4(Mx,Mxy,Qx)
    L5 = get_L5(My,Mxy,Qy)

    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            #print("Now = ",(j - 1) * N_x + i - 1)
            Es[(j - 1) * N_x + i - 1] = L3 * Get_w_sin_x(i) * Get_w_sin_y(j)
            Es[(j - 1) * N_x + i - 1 + N] = L2 * Get_w_sin_x(i) * Get_w_sin_y_2(j) #V
            Es[(j - 1) * N_x + i - 1 + N * 2] = L1 * Get_w_sin_x_2(i) * Get_w_sin_y(j) #U
            Es[(j - 1) * N_x + i - 1 + N * 3] = L4 * Get_w_cos_x(i) * Get_w_sin_y(j)  #PsiX
            Es[(j - 1) * N_x + i - 1 + N * 4] = L5 * Get_w_sin_x(i) * Get_w_cos_y(j)  #PsiY

    #print("ES",Es)
    return Es
def Get_Jacobian_Bubnov(Function_E, Result_w):
    #print("Jacodi start)")
    Jacobian = [0] * N * 5
    Def_Function = Function_E
    Symbol_Function = [0] * N * 5

    #print("Function E",Function_E)
    for i in range(0, N):
        Symbol_Function[i] = 'w' + str(i + 1)
        Symbol_Function[i + N] = 'u' + str(i + 1)
        Symbol_Function[i + N * 2] = 'v' + str(i + 1)
        Symbol_Function[i + N * 3] = 'PsiX' + str(i + 1)
        Symbol_Function[i + N * 4] = 'PsiY' + str(i + 1)

    for i in range(0, N * 5):
        Jacobian[i] = [0] * N * 5
        Def_Function[i] = Function_E[i]
        #print("D[",i,"] =",Def_Function[i])

    for row in range(0, N * 5):
        for column in range(0, N * 5):
            Jacobian[row][column] = Def_Function[row].diff(Symbol_Function[column])
            for W_coefs in range(0, N * 5):
                Jacobian[row][column] = Jacobian[row][column].subs((Symbol_Function[W_coefs]), (Result_w[W_coefs]))
            #print("J[", row, "][", column, "] =", Jacobian[row][column])
    #print("Result J =", Jacobian)
    return Jacobian
def Get_New_iterarion_Bubnov(Function_E, Jackobi_inv, Last_step_X, a):
    Def_Function = [0] * N * 5
    Symbol_Function = [0] * N * 5
    #print( "A = ", a)
    for i in range(0, N):
        Symbol_Function[i] = 'w' + str(i + 1)
        Symbol_Function[i + N] = 'u' + str(i + 1)
        Symbol_Function[i + N * 2] = 'v' + str(i + 1)
        Symbol_Function[i + N * 3] = 'PsiX' + str(i + 1)
        Symbol_Function[i + N * 4] = 'PsiY' + str(i + 1)

    for i in range(0, N * 5):
        Def_Function[i] = Function_E[i]
        #print(i, " Func e =",Function_E[i])
    #for i in range(0, N * 5):
        #for j in range(0, N * 5):
            #Def_Function[i] = Function_E[i].diff(Symbol_Function[j])
            #print(Jackobi_inv[i][j])
    for i in range(0, N * 5):
        for W_coefs in range(0, N * 5):
            Def_Function[i] = Def_Function[i].subs(Symbol_Function[W_coefs], (Last_step_X[W_coefs]))
        #print("Def_Function[", i, "] = ", Def_Function[i])

    Def_Function = sym.Matrix(Def_Function)
    W = sym.Matrix(Last_step_X)
    #print(Last_step_X)
    #print("w = ", W)
    # print("W =" ,W)
    # print("a =",a)
    # print("Ja =",Jackobi_inv)
    # print("DF =",Def_Function)
    Result = W - a * (Jackobi_inv * Def_Function)
    #W = W * Jackobi_inv
    #print("Result = ", Result)
    return Result
def Bubnov_Iter(Function_E, eps, w0):
    All_Results = []
    Res_now = w0
    Now_eps = 1
    Check_Loop = 1
    Count_Iterions = 0
    New_Eps = [0] * N*5

    while (Now_eps > eps):
        Count_Iterions += 1
        Max_eps = 0

        Jacobi = Get_Jacobian_Bubnov(Function_E, Res_now)
        print("Jacobi = ", Jacobi[0])
        Jackobi_Matrix = sym.Matrix(Jacobi)
        #print("j =",Jackobi_Matrix)
        #for i in range(5):
            #for j in range(5):
                #print("Jacobi = ", Jackobi_Matrix[i][j])

        Jacobi_Invariant = Jackobi_Matrix.inv()
        print("Jacobi = ", Jacobi_Invariant)


        Res_new = Get_New_iterarion_Bubnov(Function_E, Jacobi_Invariant, Res_now, Check_Loop)
        print("Res_now = ",Res_new)
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
        #print("R_now",Res_now)

    return Res_now
def Bubnov_Loop(Nx,Ny,Nxy,Mx,My,Mxy,Qx,Qy,Tetta_1,Tetta_2,Q_find):
    print("Start loop")
    # Q_now = 0.01
    # Q_now = 0
    q_T = 3.14 / 10
    Q_step = q_T / 2
    Miz = 0
    Q_now = 0.01
    q_for_graph = []
    w_for_graph = []
    q_for_graph.append(0)
    w_for_graph.append(0)
    W_Result = [0] * (1000)
    W_Result[0] = [0.0001] * N * 5

    Es = [0] * N * 5
    Es = Get_Es_Bub(Nx,Ny,Nxy,Mx,My,Mxy,Qx,Qy,Tetta_1,Tetta_2,Q_find)

    queue = multiprocessing.Queue()
    timer = time.time()
    p = [0] * (N * 5-N+1)

    Simp = 1
    if Simp == 1:
        for i in range(0, N * 5):
            print("Es ",i,"= ",Es[i])
            p[i] = multiprocessing.Process(target=intagrete_Es_Simp_Bubnov, args=(Es[i], queue, i,))
    if Simp == 0:
        for i in range(0, N * 5):
            p[i] = multiprocessing.Process(target=intagrete_Es, args=(Es[i], queue, i,))
    for i in range(0, N * 5):
        p[i].start()

    t = 0
    if N == 1:
        t = 100
    if N == 4:
        t = 360

    p[0].join(timeout=t)
    print("p", 0, " - join")
    for i in range(1, N * 5):
        p[i].join(timeout=t)
        print("p", i, " - join")

    j = N
    Buf_Function = [0] * N * 5
    while queue.qsize() > 0:
        buf = queue.get()
        print(buf)
        function = buf[0]
        index = buf[1]
        Buf_Function[index] = (1/2) * function
        print("j =", j, " =  ", buf)
        j += 1
    print('Time %.6f' % (time.time() - timer))

    j = 1
    while (Q_now < 6):
        Dop_L3 = Q_now
        Dop_Integrals = [0]*N
        for jj in range(1, N_x + 1):
            for i in range(1, N_y + 1):
                Dop_Integrals[(jj - 1) * N_x + i - 1] = Dop_L3 * Get_w_sin_x(i) * Get_w_sin_y(jj)  # W

        print("Q_now = ", Q_now)
        print("Start integrate 7")
        for i in range(0,N):
            Dop_Integrals[i] = integrate(Dop_Integrals[i], (Xx, Start_integral, A_lenght_x))
            Buf_Function[i] += integrate(Dop_Integrals[i], (Yy, Start_integral, B_lenght_y))
        print("End 7 ")

        W_Result[j] = Bubnov_Iter(Buf_Function, eps, W_Result[j-1])

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

        for i in range(1, N + 1):
            print("W(",j,")(",i,") = ",W_val[i - 1])
        for i in range(1, N + 1):
            print("U(", j, ")(", i, ") = ", U_val[i - 1])
        for i in range(1, N + 1):
            print("V(", j, ")(", i, ") = ", V_val[i - 1])
        for i in range(1, N + 1):
            print("PsiX(", j, ")(", i, ") = ", PsiX_val[i - 1])
        for i in range(1, N + 1):
            print("PsiY(", j, ")(", i, ") = ", PsiY_val[i - 1])

        #Miz = Draw_3d_Sigmas_main('Sigma_i', W_val, 1, U_function, V_function, W_Function, PsiX_function, PsiY_function,
                                  #z,W_val, U_val, V_val, PsiX_val, PsiY_val,Sigma_x,Sigma_y,Sigma_tay)

        New_w = Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3)
        Last_w = w_for_graph[len(w_for_graph)-1]

        #New_w = Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3)
        #print("W(l/2) = ", New_w)
        #New_w = Get_W_Plane(A_lenght_x / 4, B_lenght_y / 4, 0, W_val, 3)
        #print("W(l/4) = ", New_w)

        Q_now += Q_step
        j += 1
        print("Miz = ", Miz)

        if New_w > Last_w:
            w_for_graph.append(Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3))
            q_for_graph.append(3 * (Q_now - Q_step))
        else:
            print("trouble")
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
    print("Get_Bub")
def Bubnov_byPoint(Nx,Ny,Nxy,Mx,My,Mxy,Qx,Qy,Tetta_1,Tetta_2,Q_find):
    print("Start loop")
    # Q_now = 0.01
    # Q_now = 0
    q_T = 3.14 / 10
    Q_step = q_T / 2
    Miz = 0
    Q_now = 0.01
    q_for_graph = []
    w_for_graph = []
    q_for_graph.append(0)
    w_for_graph.append(0)
    W_Result = [0] * (1000)
    W_Result[0] = [0.0001] * N * 5

    Es = [0] * N * 5
    Es = Get_Es_Bub(Nx,Ny,Nxy,Mx,My,Mxy,Qx,Qy,Tetta_1,Tetta_2,Q_find)

    queue = multiprocessing.Queue()
    timer = time.time()
    p = [0] * (N * 5)
    Buf_Function = [0] * N * 5

    for index_thread in range (N):
        Simp = 1
        if Simp == 1:
            for i in range(index_thread * 5,index_thread * 5 + 5):
                print("Es ",i,"= ",Es[i])
                p[i] = multiprocessing.Process(target=intagrete_Es_Simp_Bubnov, args=(Es[i], queue, i,))
        if Simp == 0:
            for i in range(index_thread * 5,index_thread * 5 + 5):
                p[i] = multiprocessing.Process(target=intagrete_Es, args=(Es[i], queue, i,))
        for i in range(index_thread * 5,index_thread * 5 + 5):
            print("p ",i," start")
            p[i].start()

        t = 0
        if index_thread == 0:
            if N == 1:
                t = 100
            if N == 4:
                t = 550
        else:
            if N == 1:
                t = 100
            if N == 4:
                t = 260

        p[index_thread * 5].join()
        print("p", 0, " - join")
        for i in range(index_thread * 5+1,index_thread * 5 + 5):
            p[i].join()
            print("p", i, " - join")

        j = N
        while queue.qsize() > 0:
            buf = queue.get()
            function = buf[0]
            index = buf[1]
            Buf_Function[index] = function
            print("j =", j, " =  ", buf)
            j += 1
    print('Time %.6f' % (time.time() - timer))

    j = 1
    #Dop_L3 = Q_find
    Dop_L3 = 3.14
    Dop_Integrals = [0]*N
    for jj in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Dop_Integrals[(jj - 1) * N_x + i - 1] = Dop_L3 * Get_w_sin_x(i) * Get_w_sin_y(jj)  # W

    print("Q_now = ", Q_now)
    print("Start integrate 7")
    for i in range(0,N):
        Dop_Integrals[i] = integrate(Dop_Integrals[i], (Xx, Start_integral, A_lenght_x))
        Buf_Function[i] += integrate(Dop_Integrals[i], (Yy, Start_integral, B_lenght_y))
    print("End 7 ")

    W_Result[j] = Bubnov_Iter(Buf_Function, eps, W_Result[j-1])

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

    for i in range(1, N + 1):
        print("W(",j,")(",i,") = ",W_val[i - 1])
    for i in range(1, N + 1):
        print("U(", j, ")(", i, ") = ", U_val[i - 1])
    for i in range(1, N + 1):
        print("V(", j, ")(", i, ") = ", V_val[i - 1])
    for i in range(1, N + 1):
        print("PsiX(", j, ")(", i, ") = ", PsiX_val[i - 1])
    for i in range(1, N + 1):
        print("PsiY(", j, ")(", i, ") = ", PsiY_val[i - 1])

    #Miz = Draw_3d_Sigmas_main('Sigma_i', W_val, 1, U_function, V_function, W_Function, PsiX_function, PsiY_function,
                              #z,W_val, U_val, V_val, PsiX_val, PsiY_val,Sigma_x,Sigma_y,Sigma_tay)

    Res_w = Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3)
    print("W(l/2) =", Res_w)
    Res_w = Get_W_Plane(A_lenght_x / 4, B_lenght_y / 4, 0, W_val, 3)
    print("W(l/4) =", Res_w)

    #New_w = Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3)
    #print("W(l/2) = ", New_w)
    #New_w = Get_W_Plane(A_lenght_x / 4, B_lenght_y / 4, 0, W_val, 3)
    #print("W(l/4) = ", New_w)

    j += 1
    print("Miz = ", Miz)

    print("Get_Bub")
#Nuton
def Get_Jacobian_Nuton(Function_E):
    Jacobian = []
    Def_Function = [0] * N * 5
    Symbol_Function = [0] * N * 5

    timer = time.time()
    for i in range(0, N):
        Symbol_Function[i] = 'w' + str(i + 1)
        Symbol_Function[i + N] = 'u' + str(i + 1)
        Symbol_Function[i + N * 2] = 'v' + str(i + 1)
        Symbol_Function[i + N * 3] = 'PsiX' + str(i + 1)
        Symbol_Function[i + N * 4] = 'PsiY' + str(i + 1)
    if Multi == 1:
        Data = []
        for i in range(0, N * 5):
            Data.append((Function_E,Symbol_Function[i], i))
        with mp.Pool(N * 5) as p:
            p_result = p.map(get_diff_2, Data)
        for i in range(0, N * 5):
            Def_Function[i] = p_result[i]
        #print('First_Jacobi_multi = Time %.6f' % (time.time() - timer))
    else:
        for i in range(0, N * 5):
            Def_Function[i] = Function_E.diff(Symbol_Function[i])
        #print('First_Jacobi_not multi = Time %.6f' % (time.time() - timer))

    timer = time.time()
    Data = []
    if Multi == 1:
        for row in range(0, N * 5):
            Jacobian.append([])
            for column in range(0, N * 5):
                Jacobian[row].append([])
                Data.append((Def_Function[row], Symbol_Function[column],  row, column))
        with mp.Pool(60) as p:
            p_result = p.map(get_diff_2arg_2, Data)
        for row in range(0, N * 5):
            for column in range(0, N * 5):
                Jacobian[row][column] = p_result[row * (N * 5) + column][0]
        #print('Second_Jacobi_multi = Time %.6f' % (time.time() - timer))
    else:
        for row in range(0, N * 5):
            Jacobian.append([])
            for column in range(0, N * 5):
                Jacobian[row].append([])
                Jacobian[row][column] = Def_Function[row].diff(Symbol_Function[column])
        #print('Second_Jacobi_not_multi = Time %.6f' % (time.time() - timer))

    return Def_Function,Jacobian
def Get_New_iterarion_Nuton(Def_function_get, Jackobi_inv, Last_step_X, a):
    Def_Function = Def_function_get.copy()
    Symbol_Function = [0] * N * 5
    #print( "A = ", a)
    for i in range(0, N):
        Symbol_Function[i] = 'w' + str(i + 1)
        Symbol_Function[i + N] = 'u' + str(i + 1)
        Symbol_Function[i + N * 2] = 'v' + str(i + 1)
        Symbol_Function[i + N * 3] = 'PsiX' + str(i + 1)
        Symbol_Function[i + N * 4] = 'PsiY' + str(i + 1)

    Result_function = []
    Symbol_data = []
    Data = []
    for W_coefs in range(0, N * 5):
        Symbol_data.append((Symbol_Function[W_coefs], Last_step_X[W_coefs]))
    if Multi ==1 and N!=1:
        for i in range(0, N * 5):
            for W_coefs in range(0, N * 5):
                Def_Function[i] = Def_Function[i].subs(Symbol_Function[W_coefs], (Last_step_X[W_coefs]))
        for i in range(0, N * 5):
            Data.append((Def_Function[i], Symbol_data))
        with mp.Pool(60) as p:
            p_result = p.map(Sub_data, Data)
        for i in range(0, N * 5):
                Result_function.append(p_result[i])
    else:
        for i in range(0, N * 5):
            Result_function.append(0)
            for W_coefs in range(0, N * 5):
                    Result_function[i] = Def_Function[i].subs(Symbol_data)

    Result_function = sym.Matrix(Result_function)
    W = sym.Matrix(Last_step_X)

    step = a
    buf = Jackobi_inv * Result_function
    check =1
    new_step =0
    while(check):
        check = 0
        if abs((buf[0]*step)/4) > abs(W[0]) or (W[0] - step * buf[0]) < 0:
            step /= 2
            check = 1
            new_step = 1


    Result = W - step * buf
    print("A = ", step)
    if new_step == 1:
        step*=2
        step+=step/10
    return Result,step
def Nuton_Loop(Es_Get,W_Function_get,sigma_x,sigma_y,sigma_tay):
    print("Start loop")
    Q_now = Q_step

    Miz = []
    color = ['-.b','-b','-.r','-r','-.m','-m']
    q_for_graph = []
    w_for_graph = []
    q_for_graph_4 = []
    w_for_graph_4 = []
    for i in range(0,3):
        Miz.append([])
        Miz[i].append([])
        Miz[i].append(0)
        Miz[i].append(0)
        Miz[i].append(0)
        Miz[i].append(0)
        w_for_graph.append([])
        q_for_graph.append([])
        q_for_graph_4.append([])
        w_for_graph_4.append([])
        q_for_graph[i].append(0)
        w_for_graph[i].append(0)
        q_for_graph_4[i].append(0)
        w_for_graph_4[i].append(0)

    W_Result = [0] * (1000)
    W_Result[0] = [0.1] * N*5

    Q_function = q_function(Q_now,0)
    Es_Get[7] = (-2) * (Q_function * W_Function_get)

    timer = time.time()
    Es = Es_Get.copy()
    p_result = []
    Data = []
    Buf_Function = 0
    for i in range(0, Count_Es):
        Data.append((Es[i], i))
    with mp.Pool(Count_Es) as p:
        p_result = p.map(intagrete_Es_Simp_2, Data)
    for i in range(0, Count_Es):
        Buf_Function += (1 / 2) * p_result[i]
    print('End integrade = Time %.6f' % (time.time() - timer))

    for H_value in range(1,4):
        j=1
        index_g = 0
        Q_now = Q_step
        Miz_stop = 0
        Miz_now = 0
        new_sigma_x = sigma_x.subs([(Hh, h * H_value)])
        new_sigma_y = sigma_y.subs([(Hh, h * H_value)])
        new_sigma_tay = sigma_tay.subs([(Hh, h * H_value)])
        # Цикл по разным q
        while (Q_now < q_find*H_value):
            print("\nIteration =", j, " and h =",H_value)
            print("Q_now = ", Q_now,"\n")

            New_Buf_Function = Buf_Function.subs([('q',Q_now),(Hh,h*H_value)])

            #print("W_Result_send = ",W_Result[j-1])
            W_Result[j] = Nuton_Iter(New_Buf_Function, eps, W_Result[j-1])

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
            for i in range(1, N + 1):
                print("W(",j,")(",i,") = ",W_val[i - 1])
            for i in range(1, N + 1):
                print("U(", j, ")(", i, ") = ", U_val[i - 1])
            for i in range(1, N + 1):
                print("V(", j, ")(", i, ") = ", V_val[i - 1])
            for i in range(1, N + 1):
                print("PsiX(", j, ")(", i, ") = ", PsiX_val[i - 1])
            for i in range(1, N + 1):
                print("PsiY(", j, ")(", i, ") = ", PsiY_val[i - 1])

            Q_now += Q_step*H_value
            j += 1
            w_2 = Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3)
            w_4 = Get_W_Plane(A_lenght_x / 4, B_lenght_y / 4, 0, W_val, 3)

            #Miz_now,x_get,y_get = get_Mizz(W_val, U_val, V_val, PsiX_val, PsiY_val,new_sigma_x,new_sigma_y,new_sigma_tay)
            Miz_now, x_get, y_get = get_Mizz_2(W_val, U_val, V_val, PsiX_val, PsiY_val, new_sigma_x, new_sigma_y, new_sigma_tay)
            #print("Mizz =",Miz_now)
            #print("x_get =", x_get)
            #print("y_get =", y_get)

            if w_2 < w_for_graph[H_value-1][index_g]:
                Q_now -=(Q_step*H_value)/2
            else:
                index_g+=1
                w_for_graph[H_value-1].append(w_2)
                q_for_graph[H_value-1].append(Q_now - Q_step*H_value)
                w_for_graph_4[H_value-1].append(w_4)
                q_for_graph_4[H_value-1].append(Q_now - Q_step*H_value)
            print(" M =",Miz_now," and stop =",Miz_stop)
            if Miz_now > 1 and Miz_stop == 0 and j>3:
                Miz[H_value - 1][0] = Q_now - Q_step*H_value
                Miz[H_value - 1][1] = Miz_now
                Miz[H_value - 1][2] = x_get
                Miz[H_value - 1][3] = y_get
                print("STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                Miz_stop = 1
        if Miz_stop !=1:
            Miz[H_value - 1][0] = Q_now - Q_step * H_value
            Miz[H_value - 1][1] = Miz_now

        plt.plot(w_for_graph[H_value-1], q_for_graph[H_value-1],color[H_value*2-1])
        plt.plot(w_for_graph_4[H_value-1], q_for_graph_4[H_value-1],color[(H_value-1)*2])
    #Ww = [0] * 2
    #Qq = [0] * 2
    #Ww.append(w_for_graph[0])
    #Ww.append(w_for_graph[len(w_for_graph)-1])
    #Qq.append(q_for_graph[0])
    #Qq.append(q_for_graph[len(q_for_graph) - 1])
    for i in range(0,3):
        print("When H =",h*(i + 1 )," Miz = ",Miz[i][1]," when q =",Miz[i][0]," in point x,y =(",Miz[i][2],",",Miz[i][3],") !")
    if Change == 1:
        plt.legend(['T300/976 a in l/2', 'T300/976 a in l/4', 'T300/976 2*a in l/2', 'T300/976 a in l/4',
                    'T300/976 3*a in l/2', 'T300/976 3*a in l/4'])
    if Change == 2:
        plt.legend(['Org_glass a in l/2', 'Org_glass a in l/4', 'Org_glass 2*a in l/2', 'Org_glass a in l/4',
                    'Org_glass 3*a in l/2', 'Org_glass 3*a in l/4'])
    if Change == 3:
        plt.legend(['Still a in l/2', 'Still a in l/4', 'Still 2*a in l/2', 'Still a in l/4', 'Still 3*a in l/2',
                    'Still 3*a in l/4'])

    plt.title("График зависимости W-q")
    plt.xlabel("W [мм]")
    plt.ylabel("q [кПа]")
    plt.show()
    return W_Result
def Nuton_byPoint(Es_Get,W_Function_get):
    print("Start loop")
    print("Q_find = ", q_find)
    Q_function = q_function(q_find,0)
    Q_function = Q_function.subs('q',q_find)
    Es_Get[7] = (-2) * (Q_function * W_Function_get)
    Es = Es_Get.copy()

    W_Result = [0] * (5)
    W_Result[0] = [0.0001] * N*5
    #W_Result[0][0] = Start_w0

    timer = time.time()
    p_result =[]
    Data = []
    Buf_Function = 0
    for i in range(0,Count_Es):
        Data.append((Es[i],i))
    with mp.Pool(Count_Es) as p:
        p_result = p.map(intagrete_Es_Simp_2,Data)
    for i in range(0,Count_Es):
        Buf_Function += (1 / 2) * p_result[i]
    print('End integrade = Time %.6f' % (time.time() - timer))

    Buf_Function = Buf_Function.subs(Hh,h)
    W_Result[1] = Nuton_Iter(Buf_Function, eps, W_Result[0])

    W_val = [] * N
    U_val = [] * N
    V_val = [] * N
    PsiX_val = [] * N
    PsiY_val = [] * N

    for i in range(1, N + 1):
        W_val.append(W_Result[1][i - 1])
        U_val.append(W_Result[1][i - 1 + N])
        V_val.append(W_Result[1][i - 1 + 2 * N])
        PsiX_val.append(W_Result[1][i - 1 + 3 * N])
        PsiY_val.append(W_Result[1][i - 1 + 4 * N])
    j=0
    for i in range(1, N + 1):
        print("W(",1,")(",i,") = ",W_val[i - 1])
    for i in range(1, N + 1):
        print("U(", 1, ")(", i, ") = ", U_val[i - 1])
    for i in range(1, N + 1):
        print("V(", 1, ")(", i, ") = ", V_val[i - 1])
    for i in range(1, N + 1):
        print("PsiX(", 1, ")(", i, ") = ", PsiX_val[i - 1])
    for i in range(1, N + 1):
        print("PsiY(", 1, ")(", i, ") = ", PsiY_val[i - 1])

    Res_w_2 = Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, 0, W_val, 3)
    Res_w_4 = Get_W_Plane(A_lenght_x / 4, B_lenght_y / 4, 0, W_val, 3)
    print("W(l/2) =", Res_w_2)
    print("W(l/4) =", Res_w_4)
    return W_Result[1]
# Nut iter
def Nuton_Iter(Function_E, eps, w0):
    All_Results = []
    Result_now = w0

    Now_eps = 1
    Step_x = 1
    Count_Iterions = 0

    New_Eps = [0] * N*5
    Eps_all = 1

    #Create symbol massiv
    Symbol_Function = [0] * N * 5
    for i in range(0, N):
        Symbol_Function[i] = 'w' + str(i + 1)
        Symbol_Function[i + N] = 'u' + str(i + 1)
        Symbol_Function[i + N * 2] = 'v' + str(i + 1)
        Symbol_Function[i + N * 3] = 'PsiX' + str(i + 1)
        Symbol_Function[i + N * 4] = 'PsiY' + str(i + 1)
    Symbol_data = []
    for W_coefs in range(0, N * 5):
        Symbol_data.append((Symbol_Function[W_coefs], Result_now[W_coefs]))

    #Get DF and Jacobi
    timer = time.time()
    Def_Function,Jacobi = Get_Jacobian_Nuton(Function_E)

    #print('Get jacobi all time = Time %.6f' % (time.time() - timer))
    #Jacobi buf
    Jacobi_new = []
    for row in range(0, N * 5):
        Jacobi_new.append([])
        for column in range(0, N * 5):
            Jacobi_new[row].append(0)

    #print("Start nuton_iter")
    while (Now_eps > eps):
        timer = time.time()
        Count_Iterions += 1
        Jacobi_buf = Jacobi.copy()
        #print("Df(0) = ", Def_Function[0])
        #Subs jacobi
        if Multi == 1 and N != 1:
            Data = []
            for row in range(0, N * 5):
                for column in range(0, N * 5):
                    Data.append((Jacobi_buf[row][column], Symbol_data))
            with mp.Pool(60) as p:
                p_result = p.map(Sub_data, Data)
            for row in range(0, N * 5):
                for column in range(0, N * 5):
                    Jacobi_new[row][column] = p_result[row * (N * 5) + column]
            #print('Subs jacodi_multi = Time %.6f' % (time.time() - timer))
        else:
            for row in range(0, N * 5):
                for column in range(0, N * 5):
                    Jacobi_new[row][column] = Jacobi_buf[row][column].subs(Symbol_data)
            #print('Subs jacodi_not multi = Time %.6f' % (time.time() - timer))

        Jackobi_Matrix = sym.Matrix(Jacobi_new)
        Jacobi_Invariant = Jackobi_Matrix.inv()

        timer = time.time()
        Result_new,Step_x = Get_New_iterarion_Nuton(Def_Function, Jacobi_Invariant, Result_now, Step_x)
        #print('Get_New_iterarion = Time %.6f' % (time.time() - timer))

        print("Result_new = ",Result_new)
        np.array(Result_new).astype(np.float64)

        Max_eps = 0
        for i in range(0, N):
            New_Eps[i] = abs(Result_new[i] - Result_now[i])
            if (New_Eps[i] > Max_eps):
                Max_eps = New_Eps[i]

        Last_Result = Result_now
        Now_eps = Max_eps
        Result_now = Result_new

        All_Results.append(Result_now)
        if Count_Iterions > 25:
            Step_x = Step_x / 5
            Result_now = Last_Result
            if All_Results[Count_Iterions - 3][0] < All_Results[Count_Iterions - 2][0] < \
                    All_Results[Count_Iterions - 1][0]:
                if Step_x != 1:
                    Step_x *= 5
                Result_now = Result_new
            if All_Results[Count_Iterions - 3][0] > All_Results[Count_Iterions - 2][0] > \
                    All_Results[Count_Iterions - 1][0]:
                if Step_x != 1:
                    Step_x *= 5
                Result_now = Result_new
        if Count_Iterions > 50:
            Now_eps = 0
        #print("\n Iteration =", Count_Iterions, " and eps =", Now_eps, " \n")

    return Result_now

if __name__ == '__main__':
    # Main code
    timer = time.time()

    # Collect W-V func
    w_vals = Get_w_coefs()
    u_vals = Get_u_coefs()
    v_vals = Get_v_coefs()
    PsiX_vals = Get_PsiX_coefs()
    PsiY_vals = Get_PsiY_coefs()

    #print("Vals")
    #print(w_vals)
    #print(u_vals)
    #print(v_vals)
    #print(PsiX_vals)
    #print(PsiY_vals)

    W_function = Get_W_function_vals(w_vals)
    U_function = Get_U_function_vals(u_vals)
    V_function = Get_V_function_vals(v_vals)
    PsiX_function = Get_PsiX_function_vals(PsiX_vals)
    PsiY_function = Get_PsiY_function_vals(PsiY_vals)

    #print("Func")
    #print(W_function)
    #print(U_function)
    #print(V_function)
    #print(PsiX_function)
    #print(PsiY_function)

    # Заготовочка array's
    Max_W_values = [0] * 3
    Max_Sigmas_values = [0] * 3
    Es_main = [0] * Count_Es
    W_middle_values = []
    W_middle_2_values = []
    Q_values = []

    # Заготовочка for results
    Sigma_x_Function = 5 * Xx
    Sigma_y_Function = 5 * Xx
    Tay_xy_Function = 5 * Xx
    Q_function = 5 * Xx

    # Sigmas
    U_function_buf = U_function.copy()
    V_function_buf = V_function.copy()
    W_function_buf = W_function.copy()
    PsiX_function_buf = PsiX_function.copy()
    PsiY_function_buf = PsiY_function.copy()


    # Es main
    z_num = 0
    z_val = -h/2
    KSI_1 = ksi_1(PsiX_function,PsiY_function)
    KSI_2 = ksi_2(PsiX_function,PsiY_function)
    KSI_12 = ksi_12(PsiX_function,PsiY_function)
    TETTA_1 = Tetta_1(W_function,U_function)
    TETTA_2 = Tetta_2(W_function,V_function)
    F_EX = e_x(U_function,V_function,W_function,TETTA_1)
    F_EY = e_y(U_function,V_function,W_function,TETTA_2)
    F_XY = y_xy(U_function,V_function,W_function,TETTA_1, TETTA_2)
    F_XY_part =y_xy_part(U_function,V_function,W_function,TETTA_1, TETTA_2)


    F_EX_z = e_x_z(F_EX, KSI_1)
    F_EY_z = e_y_z(F_EY, KSI_2)
    F_XY_z = y_xy_z(F_XY, KSI_12)

    #z_val = 0
    type_z = 1
    if Change == 1:
        Sigma_x_Function = Get_Sigma_x_Orto(U_function_buf, V_function_buf, W_function_buf, E1, nu_12, nu_21,z_val, PsiX_function_buf, PsiY_function_buf)
        Sigma_y_Function = Get_Sigma_y_Orto(U_function_buf, V_function_buf, W_function_buf, E2, nu_12, nu_21, z_val, PsiX_function_buf, PsiY_function_buf)
        Tay_xy_Function = Get_Sigma_tay_Orto(U_function_buf, V_function_buf, W_function_buf, G12, z_val, PsiX_function_buf, PsiY_function_buf)
    if Change == 2 or Change == 3:
        if type_z == 1:
            Sigma_x_Function = Get_Sigma_x_Izo(E1, nu_12, z_val, F_EX_z, F_EY_z, KSI_1, KSI_2)
            Sigma_y_Function = Get_Sigma_y_Izo(E1, nu_12, z_val, F_EX_z, F_EY_z, KSI_1, KSI_2)
            Tay_xy_Function = Get_Sigma_tay_Izo(G12, nu_12, z_val, F_XY_z, KSI_12)
        if type_z == 0:
            Sigma_x_Function = Get_Sigma_x_Izo(E1, nu_12, z_val, F_EX, F_EY, KSI_1, KSI_2)
            Sigma_y_Function = Get_Sigma_y_Izo(E1, nu_12, z_val, F_EX, F_EY, KSI_1, KSI_2)
            Tay_xy_Function = Get_Sigma_tay_Izo(G12, nu_12, z_val, F_XY, KSI_12)

    #print("Sigma_x_Function = ", Sigma_x_Function)
    #print("Sigma_y_Function = ", Sigma_y_Function)
    #print("Tay_xy_Function = ", Tay_xy_Function)
    integral_type = 1
    part = 1
    if integral_type == 1:
        if Change == 1:
            Es_main[0] = N_x_Orto(z_num, U_function, V_function, W_function, E1, nu_12, nu_21) * e_x(U_function, V_function, W_function)
            Es_main[1] = N_y_Orto(z_num, U_function, V_function, W_function, E1, nu_12, nu_21) * e_y(U_function, V_function, W_function)
            Es_main[2] = N_xy_Orto(z_num, U_function, V_function, W_function, G12) * y_xy(U_function, V_function,W_function)
            Es_main[3] = M_x_Orto(W_function, E1, nu_12, nu_21) * ksi_1(PsiX_function,PsiY_function) \
                         + M_y_Orto(W_function,E2,nu_12,nu_21) * ksi_2(PsiX_function, PsiY_function)
            Es_main[4] = 2 * M_xy_Orto(W_function, G12) * ksi_12(PsiX_function, PsiY_function)
            Es_main[5] = Q_x(PsiX_function, PsiY_function, G13, W_function, U_function)
            Es_main[6] = Q_y(PsiX_function, PsiY_function, G23, W_function, U_function)
        if Change == 2 or Change == 3:
            if part ==1:
                NX_I_part = N_x_Izo_part(E1, nu_12, nu_21, F_EX, F_EY)
                NY_I_part = N_y_Izo_part(E1, nu_12, nu_21, F_EX, F_EY)
                #NXY_I = N_xy_Izo(G12, nu_12,F_XY)
                NXY_I_part = N_xy_Izo(G12, nu_12, F_XY_part)

                MX_I = M_x_Izo(E1, nu_12,KSI_1,KSI_2)
                MY_I = M_y_Izo(E1,nu_12,KSI_1,KSI_2)
                MXY_I = M_xy_Izo(G12, nu_12,KSI_12)

                Q_X = Q_x(PsiX_function,G13, TETTA_1)
                Q_Y = Q_y(PsiY_function, G23, TETTA_2)

                Es_main[0] = NX_I_part[0] * F_EX
                #print("Es_main[0] = ",Es_main[0])
                Es_main[1] = NX_I_part[1] * F_EX
                #print("Es_main[1] = ", Es_main[1])
                Es_main[2] = NY_I_part[0] * F_EY
                #print("Es_main[2] = ", Es_main[2])
                Es_main[3] = NY_I_part[1] * F_EY
                #print("Es_main[3] = ", Es_main[3])
                Es_main[4] = NXY_I_part[0]*F_XY
                #print("Es_main[4] = ", Es_main[2])
                Es_main[5] = MX_I * KSI_1 + MY_I * KSI_2 + (MXY_I + MXY_I) * KSI_12
                #print("Es_main[5] = ", Es_main[5])
                Es_main[6] = Q_X * (PsiX_function - TETTA_1) + Q_Y * (PsiY_function - TETTA_2)
                #print("Es_main[6] = ", Es_main[6])
                Es_main[8] = NXY_I_part[1] * F_XY
                #print("Es_main[8] = ", Es_main[8])
                Es_main[9] = NXY_I_part[2] * F_XY_part[0]
                #print("Es_main[9] = ", Es_main[9])
                Es_main[10] = NXY_I_part[2] * F_XY_part[1]
                #print("Es_main[10] = ", Es_main[10])
                Es_main[11] = NXY_I_part[2] * F_XY_part[2]
                #print("Es_main[11] = ", Es_main[11])
            if part == 0:
                NX_I = N_x_Izo(E1, nu_12, nu_21, F_EX, F_EY)
                NX_I_part = N_x_Izo_part(E1, nu_12, nu_21, F_EX, F_EY)
                NY_I = N_y_Izo(E1, nu_12, nu_21, F_EX, F_EY)
                NXY_I = N_xy_Izo(G12, nu_12, F_XY)

                MX_I = M_x_Izo(E1, nu_12, KSI_1, KSI_2)
                MY_I = M_y_Izo(E1, nu_12, KSI_1, KSI_2)
                MXY_I = M_xy_Izo(G12, nu_12, KSI_12)

                Q_X = Q_x(PsiX_function, G13, TETTA_1)
                Q_Y = Q_y(PsiY_function, G23, TETTA_2)

                Es_main[0] = NX_I * F_EX
                print("Es_main[0] = ", Es_main[0])
                Es_main[1] = NY_I * F_EY
                print("Es_main[1] = ", Es_main[1])
                Es_main[2] = (1 / 2) * (NXY_I + NXY_I) * F_XY
                print("Es_main[2] = ", Es_main[2])
                Es_main[3] = MX_I * KSI_1 + MY_I * KSI_2
                print("Es_main[3] = ", Es_main[3])
                Es_main[4] = (MXY_I + MXY_I) * KSI_12
                print("Es_main[4] = ", Es_main[4])
                Es_main[5] = Q_X * (PsiX_function - TETTA_1)
                print("Es_main[5] = ", Es_main[5])
                Es_main[6] = Q_Y * (PsiY_function - TETTA_2)
                print("Es_main[6] = ", Es_main[6])
                print("New_main")


        W_val = []
        U_val = []
        V_val = []
        PsiX_val = []
        PsiY_val = []
        W_values = []

        Es_main_buf = Es_main.copy()

        if One_point == 1:
            if Type_Resolve == 1:
                print("Start Nuton_point")
                W_values = Nuton_byPoint(Es_main_buf,W_function)
            if Type_Resolve == 2:
                print("Start Bubnov")
                W_values = Bubnov_byPoint(NX_I, NY_I, NXY_I, MX_I, MY_I, MXY_I, Q_X, Q_Y, TETTA_1, TETTA_2, 3.14)
        else:
            if Type_Resolve == 1:
                print("Start Nuton_loop")
                W_values = Nuton_Loop(Es_main_buf,W_function,Sigma_x_Function,Sigma_y_Function,Tay_xy_Function)
            if Type_Resolve == 2:
                print("Start Bubnov")
                W_values = Bubnov_Loop(NX_I, NY_I, NXY_I, MX_I, MY_I, MXY_I, Q_X, Q_Y, TETTA_1, TETTA_2, 3.14)

        print("W_vals")
        print(W_values)
        spend_time = time.time() - timer
        print("Prog work by = ", spend_time)

        print(W_values)
        for i in range(1, N + 1):
            W_val.append(W_values[i - 1])
            U_val.append(W_values[i - 1 + N])
            V_val.append(W_values[i - 1 + 2 * N])
            PsiX_val.append(W_values[i - 1 + 3 * N])
            PsiY_val.append(W_values[i - 1 + 4 * N])

        Symbol_Function = [0] * N * 5
        # print( "A = ", a)
        for i in range(0, N):
            Symbol_Function[i] = 'w' + str(i + 1)
            Symbol_Function[i + N] = 'u' + str(i + 1)
            Symbol_Function[i + N * 2] = 'v' + str(i + 1)
        for index in range(0, N):
            W_function = W_function.subs(Symbol_Function[index], W_values[index])
            U_function = U_function.subs(Symbol_Function[index + N], W_values[index + N])
            V_function = V_function.subs(Symbol_Function[index + N * 2], W_values[index + N * 2])
        #print("Results =")
        #print(W_function)
        #print(U_function)
        #print(V_function)
        Draw_3d_W(W_function,1,1)
        Draw_deformed_plot(W_function, U_function, V_function)

        # Prints data :
        #print("W values = ", W_val)
        #print("U values = ", U_val)
        #print("V values = ", V_val)
        #print("Psi X values = ", V_val)
        #print("Psi Y values = ", V_val)
        print('Time %.6f' % (time.time() - timer))


