from sympy import *
import sympy as sym

import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt


#Change var : 1 for  T300/976 (izo)      ;2 for Org_Glass                  ;3 for Still
Change = 1
#DATA
N_x = 4
N_y = 4
N = N_x*N_y
Q_Start = 0.02
Q_my = Q_Start

#Symbols
Xx = Symbol('x')
Yy = Symbol('y')

#Static_DATA
h=0.00022
A_lenght_x = 0.2
B_lenght_y = 0.2
R_1_x = 5
R_2_y = 3.33
A = 1
B = 1
K_x = 1/R_1_x
K_y = 1/R_2_y
z = 0

q_T = 1.34/100

#Static_DATA for variants:

#Org:
#T300/976
E1_T300 = 1.4 * (10**5)
E2_T300 = 0.97 * (10**4)
nu_T300 = 0.29
G_12_T300 = 0.55 * (10**4)
G_13_T300 = G_12_T300
G_23_T300 = 0.33 * (10**4)
F1_min_T300 = -1599
F1_max_T300 = 1517
F2_min_T300 = -253
F2_max_T300 = 46
F12_max_T300 = 41.4
Density_T300 = 1500

#Izo:
#Org_Glass_1:
E1_OG = 0.03 * (10**5)
E2_OG = E1_OG
nu_OG = 0.35
G_12_OG = 0.012 * (10**5)
G_13_OG = G_12_T300
G_23_OG = G_12_T300
Sigma_t_OG = 75
Density_OG = 1190

#Still:
E1_Still = 2.1 * (10**5)
E2_Still = E1_OG
nu_Still = 0.3
G_12_Still = 0.807 * (10**5)
G_13_Still = G_12_T300
G_23_Still = G_12_T300
Sigma_t_Still = 300
Density_Still = 7800


#Settings to integral
Start_integral = 0

#graph points
Size = 30

def Get_W_Plane(x_val,y_val,function,Values):
    W_result = 0.
    # Hard
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            W_result+= Values[(j-1)*N_x + i-1] * sin(i * x_val * m.pi / A_lenght_x) * sin(j * y_val * m.pi / B_lenght_y)

    return W_result
def Draw_3d_W(Function,Values_Result):
    x_array = []
    y_array = []
    z_array = [0]*Size

    for i in range (0,Size):
        z_array[i] = [0]*Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0,Size):
        x_array.append(i*step_x)
        y_array.append(i*step_y)

    for j in range(0,Size):
        for i in range(0,Size):
            z_array[i][j] = (Get_W_Plane(x_array[i],y_array[j],Function,Values_Result))


    z_array = np.array(z_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel(Function + ' (x,y)')
    ax.set_title(Function + ' (x,y)')
    plt.show()
def Draw_3d_Sigmas(Function, Values_Result,Type_Sigmas,U_function,V_function,W_Function,z_val,W_val,U_val,V_val):
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

    if Change == 1:
        if Type_Sigmas == 1:
            Sigma_x_Function = (Get_Sigma_x_Orto(U_function,V_function,W_Function,E1_T300,nu_T300,nu_T300,z_val))
        if Type_Sigmas == 2:
            Sigma_y_Function = (Get_Sigma_y_Orto(U_function,V_function,W_Function,E1_T300,nu_T300,nu_T300,z_val))
        if Type_Sigmas == 3:
            Tay_xy_Function = (Get_Sigma_tay_Orto(U_function,V_function,W_Function,G_12_T300,z_val))
    if Change == 2:
        if Type_Sigmas == 1:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
        if Type_Sigmas == 2:
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
        if Type_Sigmas == 3:
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
    if Change == 3:
        if Type_Sigmas == 1:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        if Type_Sigmas == 2:
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
        if Type_Sigmas == 3:
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))



    for i in range(N+1):
        Sigma_x_Function = Sigma_x_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Sigma_y_Function = Sigma_y_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Tay_xy_Function = Tay_xy_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])

    #print(Sigma_x_Function)
    #print(Sigma_y_Function)
    #print(Tay_xy_Function)

    for j in range(0, Size):
        for i in range(0, Size):
            if Type_Sigmas == 1:
                z_array[i][j] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if Type_Sigmas == 2:
                z_array[i][j] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
            if Type_Sigmas == 3:
                z_array[i][j] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])


    z_array = np.array(z_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel(Function + ' (x,y)')
    ax.set_title(Function + ' (x,y)')
    plt.show()

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



def Get_w_sin_x(i):
    return sin(i * Xx * m.pi / A_lenght_x)
def Get_w_sin_y(j):
    return sin(j * Yy * m.pi / B_lenght_y)

def Get_U_function_vals(u_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += u_vals[(j-1)*N_x + i-1] * Get_w_sin_x(i)*Get_w_sin_y(j)
    return Result
def Get_V_function_vals(v_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += v_vals[(j-1)*N_x + i-1] * Get_w_sin_x(i) * Get_w_sin_y(j)
    return Result
def Get_W_function_vals(w_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += w_vals[(j-1)*N_x + i-1] * Get_w_sin_x(i) * Get_w_sin_y(j)
    return Result

def ksi_1(W_function):
    Result = -(W_function.diff(Xx)).diff(Xx)
    return Result
def ksi_2(W_function):
    Result = -(W_function.diff(Yy)).diff(Yy)
    return Result
def ksi_12(W_function):
    Result_1 = -(W_function.diff(Xx)).diff(Yy)
    Result_2 = -(W_function.diff(Yy)).diff(Xx)
    Result = (Result_1 + Result_2)/2
    return Result

def e_x(U_function,V_function,W_Function):
    Result = U_function.diff(Xx) - K_x * W_Function
    return Result
def e_y(U_function,V_function,W_Function):
    Result = V_function.diff(Yy) - K_y * W_Function
    return Result
def y_xy(U_function,V_function,W_Function):
    Result = V_function.diff(Xx) + U_function.diff(Yy)
    return Result

def e_xz(z,U_function,V_function,W_Function):
    Result = e_x(U_function,V_function,W_Function) + z*ksi_1(W_Function)
    return Result
def e_yz(z,U_function,V_function,W_Function):
    Result = e_y(U_function,V_function,W_Function) + z*ksi_2(W_Function)
    return Result
def y_xz(z,U_function,V_function,W_Function):
    Result = y_xy(U_function,V_function,W_Function) + 2*z*ksi_12(W_Function)
    return Result

#Ortotrop functions
def N_x_Orto(z,U_function,V_function,W_Function,E_1,nu_12,nu_21):
    Result = (E_1 / (1 -nu_12*nu_21)) * h *(e_xz(z,U_function,V_function,W_Function) + nu_21 *  e_yz(z,U_function,V_function,W_Function))
    return Result
def N_y_Orto(z, U_function, V_function, W_Function,E_2,nu_12,nu_21):
    Result = (E_2 / (1 - nu_12 * nu_21)) * h * ( e_yz(z, U_function, V_function, W_Function) + nu_12 * e_xz(z, U_function, V_function, W_Function))
    return Result
def N_xy_Orto(z, U_function, V_function, W_Function,G_12):
    Result = G_12 * h * y_xz(z,U_function,V_function,W_Function)
    return Result

def M_x_Orto(W_Function,E_1,nu_12,nu_21):
    Result = (E_1 / (1 - nu_12 * nu_21)) * ((h**3)/12) * ( ksi_2(W_Function) + nu_21 * ksi_1(W_Function))
    return Result
def M_y_Orto(W_Function,E_2,nu_12,nu_21):
    Result = (E_2 / (1 - nu_12 * nu_21)) * ((h**3)/12) * ( ksi_1(W_Function) + nu_12 * ksi_2(W_Function))
    return Result
def M_xy_Orto(W_Function,G_12):
    Result = 2 * G_12 * ((h**3)/12) *  ksi_12(W_Function)
    return Result

def Get_Sigma_x_Orto(U_function,V_function,W_Function,E_1,nu_12,nu_21,z_val):
    Result = (E_1 / (1 -nu_12*nu_21))*(e_xz(z_val,U_function,V_function,W_Function) + nu_21 *  e_yz(z_val,U_function,V_function,W_Function))
    return Result
def Get_Sigma_y_Orto(U_function,V_function,W_Function,E_2,nu_12,nu_21,z_val):
    Result = (E_2 / (1 - nu_12 * nu_21)) * (e_yz(z_val, U_function, V_function, W_Function) + nu_12 * e_xz(z_val, U_function, V_function, W_Function))
    return Result
def Get_Sigma_tay_Orto(U_function,V_function,W_Function,G_12,z_val):
    Result = G_12 * y_xz(z_val,U_function,V_function,W_Function)
    return Result

#Izotrop functions
def N_x_Izo(z,U_function,V_function,W_Function,E,nu):
    Result = (E / (1 -nu*nu)) * h *(e_xz(z,U_function,V_function,W_Function) + nu *  e_yz(z,U_function,V_function,W_Function))
    return Result
def N_y_Izo(z, U_function, V_function, W_Function,E,nu):
    Result = (E / (1 - nu * nu)) * h * ( e_yz(z, U_function, V_function, W_Function) + nu * e_xz(z, U_function, V_function, W_Function))
    return Result
def N_xy_Izo(z, U_function, V_function, W_Function,E,nu):
    Result = (E*h*y_xz(z, U_function, V_function, W_Function)) / (2  + 2*nu)
    return Result

def M_x_Izo(W_Function,E,nu):
    Result = (E / (1 - nu * nu)) * ((h**3)/12) * ( ksi_2(W_Function) + nu * ksi_1(W_Function))
    return Result
def M_y_Izo(W_Function,E,nu):
    Result = (E / (1 - nu * nu)) * ((h**3)/12) * ( ksi_1(W_Function) + nu * ksi_2(W_Function))
    return Result
def M_xy_Izo(W_Function,E,nu):
    Result = E * ((h**3)/12) *  ksi_12(W_Function)/(1 + nu)
    return Result

def Get_Sigma_x_Izo(U_function,V_function,W_Function,E,nu,z_val):
    Result = (E / (1 -nu*nu))*(e_xz(z_val,U_function,V_function,W_Function) + nu *  e_yz(z_val,U_function,V_function,W_Function))
    return Result
def Get_Sigma_y_Izo(U_function,V_function,W_Function,E,nu,z_val):
    Result = (E / (1 - nu * nu)) * (e_yz(z_val, U_function, V_function, W_Function) + nu * e_xz(z_val, U_function, V_function, W_Function))
    return Result
def Get_Sigma_tay_Izo(U_function,V_function,W_Function,E,nu,z_val):
    Result = (E / (2 + 2*nu)) * y_xz(z_val,U_function,V_function,W_Function)
    return Result

def q_function(q_0,q_sv):
    Result = 0
    A1 = 0
    A_1 = [0] * 3
    A_2 = [0] * 3

    A_1[0] = 1 - ((A_lenght_x + A1) / (A_lenght_x - A1))**2
    A_1[1] = 4*(A_lenght_x + A1) / ((A_lenght_x - A1)**2)
    A_1[2] =  -4 / ((A_lenght_x - A1)**2)

    A_2[0] = 1 - ((B_lenght_y + A1) / (B_lenght_y - A1)) ** 2
    A_2[1] = 4 * (B_lenght_y + A1) / ((B_lenght_y - A1) ** 2)
    A_2[2] = -4 / ((B_lenght_y - A1) ** 2)


    Result = q_0*(A_1[0] + A_1[1]*Xx + A_1[2]*(Xx**2)) * (A_2[0] + A_2[1]*Yy + A_2[2]*(Yy**2)) + q_sv
    return Result
def Get_Answer(Es,w_coef,u_coef,v_coef):
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

    print("Start diff)")
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
        Buf[i] = nsimplify(Buf[i], tolerance=1e-7).evalf(16)

    for solution in linsolve(Buf, Buf_Symbols):
        Result = solution

    print("End diff)")
    return Result


#Main code
w_vals = Get_w_coefs()
u_vals = Get_u_coefs()
v_vals = Get_v_coefs()
W_Function = Get_W_function_vals(w_vals)
U_function = Get_U_function_vals(u_vals)
V_function = Get_V_function_vals(v_vals)

print("start claim")
Es_main = [0]*6
z = 0



if Change == 1:
    Es_main[0] = N_x_Orto(z,U_function,V_function,W_Function,E1_T300,nu_T300,nu_T300) * e_xz(z,U_function,V_function,W_Function)
    Es_main[1] = N_y_Orto(z,U_function,V_function,W_Function,E2_T300,nu_T300,nu_T300) * e_yz(z,U_function,V_function,W_Function)
    Es_main[2] = N_xy_Orto(z, U_function, V_function, W_Function,G_12_T300) * y_xz(z,U_function,V_function,W_Function)
    Es_main[3] = M_x_Orto(W_Function,E1_T300,nu_T300,nu_T300) * ksi_1(W_Function) + M_y_Orto(W_Function,E2_T300,nu_T300,nu_T300) * ksi_2(W_Function)
    Es_main[4] = 2 * M_xy_Orto(W_Function,G_12_T300) * ksi_12(W_Function)
    Es_main[5] = -2*q_function(Q_Start,Q_my)*W_Function
if Change == 2:
    Es_main[0] = N_x_Izo(z,U_function,V_function,W_Function,E1_OG,nu_OG) * e_xz(z,U_function,V_function,W_Function)
    Es_main[1] = N_y_Izo(z,U_function,V_function,W_Function,E1_OG,nu_OG) * e_yz(z,U_function,V_function,W_Function)
    Es_main[2] = (N_xy_Izo(z, U_function, V_function, W_Function,E1_Still,nu_Still) + N_xy_Izo(z, U_function, V_function, W_Function,E1_Still,nu_Still)) * y_xz(z,U_function,V_function,W_Function)/2
    Es_main[3] = M_x_Izo(W_Function,E1_OG,nu_OG) * ksi_1(W_Function) + M_y_Izo(W_Function,E1_OG,nu_OG) * ksi_2(W_Function)
    Es_main[4] = (M_xy_Izo(W_Function,E1_OG,nu_OG) + M_xy_Izo(W_Function,E1_OG,nu_OG)) * ksi_12(W_Function)
    Es_main[5] = -2*q_function(Q_Start,Q_my)*W_Function
if Change == 3:
    Es_main[0] = N_x_Izo(z, U_function, V_function, W_Function,E1_Still,nu_Still) * e_xz(z, U_function, V_function, W_Function)
    Es_main[1] = N_y_Izo(z, U_function, V_function, W_Function,E1_Still,nu_Still) * e_yz(z, U_function, V_function, W_Function)
    Es_main[2] = (N_xy_Izo(z, U_function, V_function, W_Function,E1_Still,nu_Still) + N_xy_Izo(z, U_function, V_function, W_Function,E1_Still,nu_Still)) * y_xz(z,U_function,V_function,W_Function)/2
    Es_main[3] = M_x_Izo(W_Function,E1_Still,nu_Still) * ksi_1(W_Function) + M_y_Izo(W_Function,E1_Still,nu_Still) * ksi_2(W_Function)
    Es_main[4] = (M_xy_Izo(W_Function,E1_Still,nu_Still) + M_xy_Izo(W_Function,E1_Still,nu_Still)) * ksi_12(W_Function)
    Es_main[5] = -2*q_function(Q_Start,Q_my) * W_Function


z= h/2
print("Start get result)")
W_values = Get_Answer(Es_main, w_vals, u_vals, v_vals)

W_val = []
U_val = []
V_val = []

print(W_values)
for i in range(1,N + 1):
    W_val.append(W_values[i - 1])
    U_val.append(W_values[i - 1 + N])
    V_val.append(W_values[i - 1 + 2*N])

#Prints data :
print(W_val)
print(U_val)
print(V_val)


#Graph's :
Draw_3d_W('W',W_val)
Draw_3d_W('U',U_val)
Draw_3d_W('V',V_val)

z=h/2

Draw_3d_Sigmas('Sigma_x',W_val,1,U_function,V_function,W_Function,z,W_val,U_val,V_val)
Draw_3d_Sigmas('Sigma_Y',W_val,2,U_function,V_function,W_Function,z,W_val,U_val,V_val)
Draw_3d_Sigmas('Tay_xy',W_val,3,U_function,V_function,W_Function,z,W_val,U_val,V_val)

