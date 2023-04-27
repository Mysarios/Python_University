import sympy as sym
from sympy import *
import numpy as np
import math as m
import matplotlib as mpl
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

#Dynamic Data
n_1 = 2
n_2 = 2
N = n_1*n_2
#Static DATA
R_1 = 20.25
R_2 = 20.25
kx = 1/R_1
ky = 1/R_2
E_1 = 2.1*(10**5)
E_2 = 2.1*(10**5)
q = 1.34/100
A = 1
B = 1
h = 0.12
nu_12 = 0.3
nu_21 = 0.3
G_12 = 0.807*(10**5)
#Symbols
Xx = Symbol('x')
Yy = Symbol('y')
#graph points
size_Graph = 20
Size = 30
Max_q_T = 2
Mm=10**5
eps = 0.001
nu = 0.3
#Settings to integral
a = A
b = B

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
    return sin(i * Xx * m.pi / A)
def Get_w_sin_y(j):
    return sin(j * Yy * m.pi / B)
def Get_w_cos(i):
    return 1 - cos(2 * i * Xx * m.pi/A)
def Get_w_cos(j):
    return 1 - cos(2 * j * Yy * m.pi/B)

def Get_U_function(w_vals):
    Result = 0
    for j in range(1, n_1 + 1):
        for i in range(1, n_2 + 1):
            Result += w_vals[(j-1)*n_1 + i-1] * Get_w_sin_x(i)*Get_w_sin_y(j)
    return Result
def Get_V_function(w_vals):
    Result = 0
    for j in range(1, n_1 + 1):
        for i in range(1, n_2 + 1):
            Result += w_vals[(j-1)*n_1 + i-1] * Get_w_sin_x(i) * Get_w_sin_y(j)
    return Result
def Get_W_function(w_vals):
    Result = 0
    for j in range(1, n_1 + 1):
        for i in range(1, n_2 + 1):
            Result += w_vals[(j-1)*n_1 + i-1] * Get_w_sin_x(i) * Get_w_sin_y(j)
    return Result

def ksi_1(W_function):
    Result = ( -1/ A ) * (W_function.diff(Xx)).diff(Xx)# - (1 / (A*B)) * ((W_function.diff(Yy)/B ) * A.diff(Yy))
    return Result
def ksi_2(W_function):
    Result = (-1 / B) * (W_function.diff(Yy)).diff(Yy)# - (1 / (A * B)) * ((W_function.diff(Xx) / A) * B.diff(Xx))
    return Result
def ksi_12(W_function):
    Result_1 = ( -1/ A ) * (W_function.diff(Xx)).diff(Yy)# - (1 / (A*B)) * ((W_function.diff(Yy)/B ) * A.diff(Yy))
    Result_2 = (-1 / B) * (W_function.diff(Xx)).diff(Yy)
    Result = Result_1 + Result_2# + (1 / (A * B)) * ((W_function.diff(Xx) / A) * A.diff(Yy) + (W_function.diff(Yy) / B) * B.diff(Xx))
    return Result

def e_x(U_function,V_function,W_Function):
    Result = ( 1/A ) * U_function.diff(Xx) - kx * W_Function # + 1 / (A*B) * A.diff(Yy) * V_function
    return Result
def e_y(U_function,V_function,W_Function):
    Result = (1 / B) * V_function.diff(Yy) - ky * W_Function #+ 1 / (A * B) * B.diff(Xx) * U_function
    return Result
def y_xy(U_function,V_function,W_Function):
    Result = (1 / A) * V_function.diff(Xx) +(1 / B) * U_function.diff(Yy)# - 1 / (A * B) *( A.diff(Yy) * U_function - B.diff(Xx) * V_function  )
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

def N_x(z,U_function,V_function,W_Function):
    Result = E_1 / (1 -nu_12*nu_21) * h *(e_xz(z,U_function,V_function,W_Function) + nu_21 *  e_yz(z,U_function,V_function,W_Function))
    return Result
def N_y(z, U_function, V_function, W_Function):
    Result = E_2 / (1 - nu_12 * nu_21) * h * ( e_yz(z, U_function, V_function, W_Function) + nu_12 * e_xz(z, U_function, V_function, W_Function))
    return Result
def N_xy(z, U_function, V_function, W_Function):
    Result = G_12 * h * y_xz(z,U_function,V_function,W_Function)
    return Result

def M_x(W_Function):
    Result = E_1 / (1 - nu_12 * nu_21) * ((h**3)/12) * ( ksi_2(W_Function) + nu_21 * ksi_1(W_Function))
    return Result
def M_y(W_Function):
    Result = E_2 / (1 - nu_12 * nu_21) * ((h**3)/12) * ( ksi_1(W_Function) + nu_12 * ksi_2(W_Function))
    return Result
def M_xy(W_Function):
    Result = 2 * G_12 * ((h**3)/12) *  ksi_12(W_Function)
    return Result

def Get_Answer(Es,w_coef,u_coef,v_coef):
    Es[0] = integrate(Es[0], (Xx, 0, a))
    Es[1] = integrate(Es[1], (Xx, 0, a))
    Es[2] = integrate(Es[2], (Xx, 0, a))
    Es[3] = integrate(Es[3], (Xx, 0, a))
    Es[4] = integrate(Es[4], (Xx, 0, a))
    Es[5] = integrate(Es[5], (Xx, 0, a))

    Es[0] = (1/2) * integrate(Es[0], (Yy, 0, b))
    Es[1] = (1/2) * integrate(Es[1], (Yy, 0, b))
    Es[2] = (1/2) * integrate(Es[2], (Yy, 0, b))
    Es[3] = (1/2) * integrate(Es[3], (Yy, 0, b))
    Es[4] = (1/2) * integrate(Es[4], (Yy, 0, b))
    Es[5] = (1/2) * integrate(Es[5], (Yy, 0, b))

    Result = []
    Result_buf = []
    Buf = []
    Zeroes = [0]*N*3
    Buf_Symbols = []

    print("Start diff)")
    for index in range(1,4):
        Result.append(Result_buf)

    for i in range(1, N + 1):
        Buf.append(Es[0].diff(w_coef[i - 1]) + Es[1].diff(w_coef[i - 1]) + Es[2].diff(w_coef[i - 1]) + Es[3].diff(w_coef[i - 1]) + Es[4].diff(w_coef[i - 1]) + Es[5].diff(w_coef[i - 1]))
        Buf_Symbols.append(w_coef[i - 1])


    for i in range(1, N + 1):
        Buf.append(Es[0].diff(u_coef[i - 1]) + Es[1].diff(u_coef[i - 1]) + Es[2].diff(u_coef[i - 1]) + Es[3].diff(w_coef[i - 1]) + Es[4].diff(u_coef[i - 1]) + Es[5].diff(u_coef[i - 1]))
        Buf_Symbols.append(u_coef[i - 1])


    for i in range(1, N + 1):
        Buf.append(Es[0].diff(v_coef[i - 1]) + Es[1].diff(v_coef[i - 1]) + Es[2].diff(v_coef[i - 1]) + Es[3].diff(v_coef[i - 1]) + Es[4].diff(v_coef[i - 1]) + Es[5].diff(v_coef[i - 1]))
        Buf_Symbols.append(v_coef[i - 1])

    for i in range(len(Buf)):
        Buf[i] = nsimplify(Buf[i],tolerance=1e-5).evalf(6)


    for solution in nonlinsolve(Buf,Buf_Symbols):
        Result = solution

    return Result

w_vals = Get_w_coefs()
u_vals = Get_u_coefs()
v_vals = Get_v_coefs()
W_Function = Get_W_function(w_vals)
U_function = Get_U_function(u_vals)
V_function = Get_V_function(v_vals)
print("start claim")
z = h/2

Es_main = [0]*6
Es_main[0] = N_x(z,U_function,V_function,W_Function) * e_xz(z,U_function,V_function,W_Function)
Es_main[1] = N_y(z,U_function,V_function,W_Function) * e_yz(z,U_function,V_function,W_Function)
Es_main[2] = (1/2) * ( N_xy(z, U_function, V_function, W_Function) * N_xy(z, U_function, V_function, W_Function)) *y_xz(z,U_function,V_function,W_Function)
Es_main[3] = M_x(W_Function) * ksi_1(W_Function) + M_y(W_Function) * ksi_2(W_Function)
Es_main[4] = (M_xy(W_Function) * M_xy(W_Function)) * ksi_12(W_Function)
Es_main[5] = -q*W_Function


print("Start get)")
W_values = Get_Answer(Es_main,w_vals,u_vals,v_vals)
print(W_values)

