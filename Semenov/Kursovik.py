from sympy import *
import numpy as np
import math as m
import matplotlib.pyplot as plt


#Change var : 1 for  T300/976 (izo)  ;2 for Org_Glass  ;3 for Still
Change = 1
Only_Value_in_Break_point = 1
Q_start = 0.018036     #1
#Q_start = 0.004655 #2
#Q_start = 0.025 #3
#Q_start = 0.000
#DATA


V = 0.2*0.2*0.00022 #m^3
P1=1500
P2=7800
P3=1190
Square = 20*20


N_x = 1
N_y = 1
N = N_x*N_y
Q_max = 245/(10**6) #1 кг максимум
mass=0
if Change == 1:
    mass = P1*V
if Change == 2:
    mass = P2*V
if Change == 3:
    mass = P3*V
Q_my = (mass/Square)/(10**6)

#Symbols
Xx = Symbol('x')
Yy = Symbol('y')

#Static_DATA
h=0.00022
#h=0.00044
#h=0.00088
A_lenght_x = 0.2
B_lenght_y = 0.2
R_1_x = 5
R_2_y = 3.33
A = 1
B = 1
K_x = 1/R_1_x
K_y = 1/R_2_y
z = 0


#Static_DATA for variants:

#Ort:
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
Sigma_t_T300 = 100

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
Size = 50

def Get_W_Plane(x_val,y_val,function,Values,type):
    W_result = 0

    # Hard
    if type == 1:
        for j in range(1, N_x + 1):
            for i in range(1, N_y + 1):
                W_result += Values[(j-1)*N_x + i-1] * sin(2*i * x_val * m.pi / A_lenght_x) * sin((2*j - 1) * y_val * m.pi / B_lenght_y)
    if type == 2:
        for j in range(1, N_x + 1):
            for i in range(1, N_y + 1):
                W_result += Values[(j-1)*N_x + i-1] * sin((2*i - 1) * x_val * m.pi / A_lenght_x) * sin(2*j * y_val * m.pi / B_lenght_y)
    if type == 3:
        for j in range(1, N_x + 1):
            for i in range(1, N_y + 1):
                W_result += Values[(j-1)*N_x + i-1] * sin((2*i - 1) * x_val * m.pi / A_lenght_x) * sin((2*j - 1) * y_val * m.pi / B_lenght_y)

    #print(W_result)
    return W_result
def Draw_3d_W(Function,Values_Result,type_f):
    x_array = []
    y_array = []
    z_array = [0]*Size
    Max_Value = 0
    #print(Values_Result)

    #print()
    #print(Function)

    for i in range (0,Size):
        z_array[i] = [0]*Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0,Size):
        x_array.append(i*step_x)
        y_array.append(i*step_y)

    for j in range(0,Size):
        for i in range(0,Size):
            z_array[i][j] = (Get_W_Plane(x_array[i],y_array[j],Function,Values_Result,type_f))
            if abs(z_array[i][j]) > abs(Max_Value):
                Max_Value = z_array[i][j]

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
    return Max_Value
def Draw_3d_Q(Function,Values_Result,type_f,q,Q_function):
    x_array = []
    y_array = []
    z_array = [0]*Size
    Max_Value = 0
    #print(Values_Result)

    #print()
    #print(Function)

    for i in range (0,Size):
        z_array[i] = [0]*Size

    step_x = A_lenght_x / (Size - 1)
    step_y = B_lenght_y / (Size - 1)

    for i in range(0,Size):
        x_array.append(i*step_x)
        y_array.append(i*step_y)

    for j in range(0,Size):
        for i in range(0,Size):
            z_array[i][j] = Q_function.subs([(Xx,x_array[i]),(Yy,y_array[j])])


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
            Sigma_x_Function = (Get_Sigma_x_Orto(U_function,V_function,W_Function,E1_T300,nu_T300,nu_T300,z_val))
        if Type_Sigmas == 2:
            Sigma_y_Function = (Get_Sigma_y_Orto(U_function,V_function,W_Function,E2_T300,nu_T300,nu_T300,z_val))
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
        Tay_xy_Function  =  Tay_xy_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])

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
            if abs(z_array[i][j]) > abs(Max_Value):
                Max_Value = z_array[i][j]


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
    return Max_Value
def Draw_3d_Sigmas_main(Function, Values_Result,Type_Sigmas,U_function,V_function,W_Function,z_val,W_val,U_val,V_val):
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
            Sigma_x_Function = (Get_Sigma_x_Orto(U_function,V_function,W_Function,E1_T300,nu_T300,nu_T300,z_val))
            Sigma_y_Function = (Get_Sigma_y_Orto(U_function,V_function,W_Function,E2_T300,nu_T300,nu_T300,z_val))
            Tay_xy_Function = (Get_Sigma_tay_Orto(U_function,V_function,W_Function,G_12_T300,z_val))
    if Change == 2:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
    if Change == 3:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))


    for i in range(N+1):
        Sigma_x_Function = Sigma_x_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Sigma_y_Function = Sigma_y_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Tay_xy_Function  =  Tay_xy_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])

    #print(Sigma_x_Function)
    #print(Sigma_y_Function)
    #print(Tay_xy_Function)


    for j in range(0, Size):
        for i in range(0, Size):
                Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                z_array[i][j] =((Max_Sigmas_values[0] ** 2) + (Max_Sigmas_values[1] ** 2) - Max_Sigmas_values[0] * Max_Sigmas_values[1] + 3 * (Max_Sigmas_values[2] ** 2)) ** (1 / 2)



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
    return Max_Value
def Draw_3d_Sigmas_main_Orto(Function, Values_Result,Type_Sigmas,U_function,V_function,W_Function,z_val,W_val,U_val,V_val):
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
            Sigma_x_Function = (Get_Sigma_x_Orto(U_function,V_function,W_Function,E1_T300,nu_T300,nu_T300,z_val))
            Sigma_y_Function = (Get_Sigma_y_Orto(U_function,V_function,W_Function,E2_T300,nu_T300,nu_T300,z_val))
            Tay_xy_Function = (Get_Sigma_tay_Orto(U_function,V_function,W_Function,G_12_T300,z_val))
    if Change == 2:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
    if Change == 3:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))


    for i in range(N+1):
        Sigma_x_Function = Sigma_x_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Sigma_y_Function = Sigma_y_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Tay_xy_Function  =  Tay_xy_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])

    #print(Sigma_x_Function)
    #print(Sigma_y_Function)
    #print(Tay_xy_Function)
    Buf_1 = (1 / F1_max_T300) + (1 / F1_min_T300)
    Buf_2 = (1 / F2_max_T300) + (1 / F2_min_T300)
    Buf_3 = 1 / (F1_max_T300 * F1_min_T300)
    Buf_4 = 1 / (F2_max_T300 * F2_min_T300)
    Buf_6 = 1 / (F12_max_T300*F12_max_T300)

    First_val = (1/2) * (Buf_4 - Buf_3)
    Second_Val = (1/2) * (Buf_3 - Buf_4)
    Third_Val = (1/2) * (Buf_3 + Buf_4)

    for j in range(0, Size):
        for i in range(0, Size):
                Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, x_array[i]), (Yy, y_array[j])])
                z_array[i][j] = First_val*(Max_Sigmas_values[0]**2) + Second_Val*(Max_Sigmas_values[1]**2) \
                                + Buf_1 * Max_Sigmas_values[0] + Buf_2 * Max_Sigmas_values[1] \
                                - Third_Val*((Max_Sigmas_values[0] - Max_Sigmas_values[2] )**2) - Buf_6*(Max_Sigmas_values[2]**2)
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
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_zlabel(Function + ' (x,y)')
        ax.set_title(Function + ' (x,y)')
        plt.show()
    return Max_Value
def Get_Sigmas_max_values(Function, Values_Result,Type_Sigmas,Sigma_x_Function,Sigma_y_Function,Tay_xy_Function,z_val,W_val,U_val,V_val):
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

    for i in range(N+1):
        if Type_Sigmas == 1:
            Sigma_x_Function = Sigma_x_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        if Type_Sigmas == 2:
            Sigma_y_Function = Sigma_y_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        if Type_Sigmas == 3:
            Tay_xy_Function  =  Tay_xy_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])

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
def Get_Sigmas_max_values_Orto(Function, Values_Result,Type_Sigmas,U_function,V_function,W_Function,z_val,W_val,U_val,V_val):
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
            Sigma_x_Function = (Get_Sigma_x_Orto(U_function,V_function,W_Function,E1_T300,nu_T300,nu_T300,z_val))
            Sigma_y_Function = (Get_Sigma_y_Orto(U_function,V_function,W_Function,E2_T300,nu_T300,nu_T300,z_val))
            Tay_xy_Function = (Get_Sigma_tay_Orto(U_function,V_function,W_Function,G_12_T300,z_val))
    if Change == 2:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function,V_function,W_Function,E1_OG,nu_OG,z_val))
    if Change == 3:
            Sigma_x_Function = (Get_Sigma_x_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
            Sigma_y_Function = (Get_Sigma_y_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))
            Tay_xy_Function = (Get_Sigma_tay_Izo(U_function, V_function, W_Function, E1_Still, nu_Still, z_val))


    for i in range(N+1):
        Sigma_x_Function = Sigma_x_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Sigma_y_Function = Sigma_y_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])
        Tay_xy_Function  =  Tay_xy_Function.subs([('w' + str(i),W_val[i-1]),('u' + str(i),U_val[i-1]),('v' + str(i),V_val[i-1])])

    #print(Sigma_x_Function)
    #print(Sigma_y_Function)
    #print(Tay_xy_Function)
    Buf_1 = (1 / F1_max_T300) + (1 / F1_min_T300)
    Buf_2 = (1 / F2_max_T300) + (1 / F2_min_T300)
    Buf_3 = 1 / (F1_max_T300 * F1_min_T300)
    Buf_4 = 1 / (F2_max_T300 * F2_min_T300)
    Buf_6 = 1 / (F12_max_T300*F12_max_T300)

    First_val = (1/2) * (Buf_4 - Buf_3)
    Second_Val = (1/2) * (Buf_3 - Buf_4)
    Third_Val = (1/2) * (Buf_3 + Buf_4)


    Max_Sigmas_values[0] = Sigma_x_Function.subs([(Xx,A_lenght_x/2), (Yy,B_lenght_y/2)])
    Max_Sigmas_values[1] = Sigma_y_Function.subs([(Xx, A_lenght_x/2), (Yy,B_lenght_y/2)])
    Max_Sigmas_values[2] = Tay_xy_Function.subs([(Xx, A_lenght_x/2), (Yy, B_lenght_y/2)])
    Max_Value = First_val*(Max_Sigmas_values[0]**2) + Second_Val*(Max_Sigmas_values[1]**2) \
                                + Buf_1 * Max_Sigmas_values[0] + Buf_2 * Max_Sigmas_values[1] \
                                - Third_Val*((Max_Sigmas_values[0] - Max_Sigmas_values[2] )**2) - Buf_6*(Max_Sigmas_values[2]**2)


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

def Get_w_sin_x(i):
    return sin((2*i-1) * Xx * m.pi / A_lenght_x)
def Get_w_sin_y(j):
    return sin((2*j-1) * Yy * m.pi / B_lenght_y)
def Get_w_sin_x_2(i):
    return sin((2*i) * Xx * m.pi / A_lenght_x)
def Get_w_sin_y_2(j):
    return sin((2*j) * Yy * m.pi / B_lenght_y)

def Get_U_function_vals(u_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += u_vals[(j-1)*N_x + i-1] * Get_w_sin_x_2(i)*Get_w_sin_y(j)
    return Result
def Get_V_function_vals(v_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += v_vals[(j-1)*N_x + i-1] * Get_w_sin_x(i) * Get_w_sin_y_2(j)
    return Result
def Get_W_function_vals(w_vals):
    Result = 0
    for j in range(1, N_x + 1):
        for i in range(1, N_y + 1):
            Result += w_vals[(j-1)*N_x + i-1] * Get_w_sin_x(i) * Get_w_sin_y(j)
    return Result

def ksi_1(W_function):
    Buf = W_function.diff(Xx)
    Result = Buf.diff(Xx) * (-1)
    return Result
def ksi_2(W_function):
    Buf = W_function.diff(Yy)
    Result = Buf.diff(Yy) * (-1)
    return Result
def ksi_12(W_function):
    Result_1 = -(W_function.diff(Xx)).diff(Yy)
    Result_2 = -(W_function.diff(Yy)).diff(Xx)
    #Result = (Result_1 + Result_2)/2
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

def e_xz(z_value,U_function,V_function,W_Function):
    Result = e_x(U_function,V_function,W_Function) + z_value*ksi_1(W_Function)
    return Result
def e_yz(z_value,U_function,V_function,W_Function):
    Result = e_y(U_function,V_function,W_Function) + z_value*ksi_2(W_Function)
    return Result
def y_xz(z_value,U_function,V_function,W_Function):
    Result = y_xy(U_function,V_function,W_Function) + 2*z_value*ksi_12(W_Function)
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
    Result = (E_1 / (1 - nu_12 * nu_21)) * (e_xz(z_val, U_function, V_function, W_Function) + nu_21 * e_yz(z_val, U_function, V_function, W_Function))
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
def Get_Answer(Es_Get,w_coef,u_coef,v_coef):
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

    #print("Start diff)")
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
    #print("End diff)")
    return Result


#Main code
#Collect W-V func
w_vals = Get_w_coefs()
u_vals = Get_u_coefs()
v_vals = Get_v_coefs()
W_Function = Get_W_function_vals(w_vals)
U_function = Get_U_function_vals(u_vals)
V_function = Get_V_function_vals(v_vals)

#Заготовочка array's
Max_W_values = [0]*3
Max_Sigmas_values = [0]*3
Es_main = [0]*6
W_middle_values = []
Q_values = []

#Заготовочка
Sigma_x_Function = 5 * Xx
Sigma_y_Function = 5 * Xx
Tay_xy_Function = 5 * Xx
Q_function = 5 * Xx


#Sigma max
QQ = 0
Check = 1
if Change == 1:
    QQ = 1
if Change == 2:
    QQ = Sigma_t_OG
if Change == 3:
    QQ = Sigma_t_Still

#Sigmas
U_function_buf = U_function.copy()
V_function_buf = V_function.copy()
W_function_buf = W_Function.copy()
z_val = -h / 2
if Change == 1:
    Sigma_x_Function = Get_Sigma_x_Orto(U_function_buf, V_function_buf, W_function_buf, E1_T300, nu_T300, nu_T300, z_val)
    Sigma_y_Function = Get_Sigma_y_Orto(U_function_buf, V_function_buf, W_function_buf, E2_T300, nu_T300, nu_T300, z_val)
    Tay_xy_Function = Get_Sigma_tay_Orto(U_function_buf, V_function_buf, W_function_buf, G_12_T300, z_val)
if Change == 2:
    Sigma_x_Function = (Get_Sigma_x_Izo(U_function_buf, V_function_buf, W_function_buf, E1_OG, nu_OG, z_val))
    Sigma_y_Function = (Get_Sigma_y_Izo(U_function_buf, V_function_buf, W_function_buf, E1_OG, nu_OG, z_val))
    Tay_xy_Function = (Get_Sigma_tay_Izo(U_function_buf, V_function_buf, W_function_buf, E1_OG, nu_OG, z_val))
if Change == 3:
    Sigma_x_Function = (Get_Sigma_x_Izo(U_function_buf, V_function_buf, W_function_buf, E1_Still, nu_Still, z_val))
    Sigma_y_Function = (Get_Sigma_y_Izo(U_function_buf, V_function_buf, W_function_buf, E1_Still, nu_Still, z_val))
    Tay_xy_Function = (Get_Sigma_tay_Izo(U_function_buf, V_function_buf, W_function_buf, E1_Still, nu_Still, z_val))

Count_num = 0
Q_now =Q_start

while Check:
    if Q_now > 0.018036:
        Check = 0
        print(W_middle_values)
        print(W_middle_values)
    Count_num +=1
    W_val = []
    U_val = []
    V_val = []
    W_values = []

    # Es main
    z_num = 0
    if Change == 1:
        Es_main[0] = N_x_Orto(z_num, U_function, V_function, W_Function, E1_T300, nu_T300, nu_T300) * e_xz(z_num, U_function, V_function,W_Function)
        Es_main[1] = N_y_Orto(z_num, U_function, V_function, W_Function, E2_T300, nu_T300, nu_T300) * e_yz(z_num, U_function,V_function,W_Function)
        Es_main[2] = N_xy_Orto(z_num, U_function, V_function, W_Function, G_12_T300) * y_xz(z_num, U_function, V_function,W_Function)
        Es_main[3] = M_x_Orto(W_Function, E1_T300, nu_T300, nu_T300) * ksi_1(W_Function) + M_y_Orto(W_Function, E2_T300,nu_T300,nu_T300) * ksi_2(W_Function)
        Es_main[4] = 2 * M_xy_Orto(W_Function, G_12_T300) * ksi_12(W_Function)
    if Change == 2:
        Es_main[0] = N_x_Izo(z_num, U_function, V_function, W_Function, E1_OG, nu_OG) * e_xz(z_num, U_function, V_function,W_Function)
        Es_main[1] = N_y_Izo(z_num, U_function, V_function, W_Function, E1_OG, nu_OG) * e_yz(z_num, U_function, V_function,W_Function)
        Es_main[2] = (N_xy_Izo(z_num, U_function, V_function, W_Function, E1_Still, nu_Still) + N_xy_Izo(z_num, U_function,V_function,W_Function, E1_Still,nu_Still)) * y_xz(z_num, U_function, V_function, W_Function)
        Es_main[3] = M_x_Izo(W_Function, E1_OG, nu_OG) * ksi_1(W_Function) + M_y_Izo(W_Function, E1_OG, nu_OG) * ksi_2(W_Function)
        Es_main[4] = (M_xy_Izo(W_Function, E1_OG, nu_OG) + M_xy_Izo(W_Function, E1_OG, nu_OG)) * ksi_12(W_Function)
    if Change == 3:
        Es_main[0] = N_x_Izo(z_num, U_function, V_function, W_Function, E1_Still, nu_Still) * e_xz(z_num, U_function,V_function, W_Function)
        Es_main[1] = N_y_Izo(z_num, U_function, V_function, W_Function, E1_Still, nu_Still) * e_yz(z_num, U_function,V_function, W_Function)
        Es_main[2] = (N_xy_Izo(z_num, U_function, V_function, W_Function, E1_Still, nu_Still) + N_xy_Izo(z_num, U_function,V_function,W_Function,E1_Still,nu_Still)) * y_xz(z_num, U_function, V_function, W_Function)
        Es_main[3] = M_x_Izo(W_Function, E1_Still, nu_Still) * ksi_1(W_Function) + M_y_Izo(W_Function, E1_Still,nu_Still) * ksi_2(W_Function)
        Es_main[4] = (M_xy_Izo(W_Function, E1_Still, nu_Still) + M_xy_Izo(W_Function, E1_Still, nu_Still)) * ksi_12(W_Function)

    Es_main[5] = (-2)*q_function(Q_now,Q_my)*W_Function
    Es_main_buf = Es_main.copy()
    Q_function = q_function(Q_now, Q_my)


    W_values = Get_Answer(Es_main_buf, w_vals, u_vals, v_vals)


    print(W_values)
    for i in range(1, N + 1):
        W_val.append(W_values[i - 1])
        U_val.append(W_values[i - 1 + N])
        V_val.append(W_values[i - 1 + 2 * N])

    for i in range(N + 1):
        Q_function = Q_function.subs('w' + str(i), W_val[i - 1])
        #Q_function = Q_function.subs(Xx,A_lenght_x/2)
        #Q_function = Q_function.subs(Yy,B_lenght_y/ 2)
    #Q_values.append(Q_function)


    z = -h/2
    Max_Sigmas_values[0] = Get_Sigmas_max_values('Sigma_x', W_val, 1, Sigma_x_Function, Sigma_y_Function, Tay_xy_Function, z, W_val, U_val, V_val)
    Max_Sigmas_values[1] = Get_Sigmas_max_values('Sigma_Y', W_val, 2, Sigma_x_Function, Sigma_y_Function, Tay_xy_Function, z, W_val, U_val, V_val)
    Max_Sigmas_values[2] = Get_Sigmas_max_values('Tay_xy',  W_val, 3, Sigma_x_Function, Sigma_y_Function, Tay_xy_Function, z, W_val, U_val, V_val)

    Sigma_max_global=0

    if Change == 1:
        Sigma_max_global = abs(Get_Sigmas_max_values_Orto('Sigma_i', W_val, 1, U_function, V_function, W_Function, z, W_val, U_val, V_val))
        #print(Sigma_max_global)
    else:
        Sigma_max_global = ((Max_Sigmas_values[0]**2) + (Max_Sigmas_values[1]**2) - Max_Sigmas_values[0]*Max_Sigmas_values[1] + 3*(Max_Sigmas_values[2]**2))**(1/2)


    W_middle_values.append(Get_W_Plane(A_lenght_x / 2, B_lenght_y / 2, W_Function, W_val,3))
    Q_now += Q_max*40

    if Sigma_max_global > QQ and Only_Value_in_Break_point != 1:
        Check = 0
    if Only_Value_in_Break_point and Count_num == 1:
        Check = 0
        Q_now -= Q_max

#Q_graph
#plt.plot(W_middle_values,Q_values)
#plt.ylabel("q , МПа")
#plt.xlabel('W(l/2,l/2) , м')
#plt.show()
print((Q_now + Q_my)*(10**6))
print(Count_num)

#Prints data :
print("W values = ",W_val)
print("U values = ",U_val)
print("V values = ",V_val)


#Graph's :

Draw_3d_Q('q ',W_val,0,Q_now,Q_function)
Max_W_values[0] = Draw_3d_W('W',W_val,3)
Max_W_values[1] = Draw_3d_W('U',U_val,1)
Max_W_values[2] = Draw_3d_W('V',V_val,2)

print("max W val = ",Max_W_values[0])
print("max U val = ",Max_W_values[1])
print("max V val = ",Max_W_values[2])


z=-h/2

if Change == 1:
    Max = Draw_3d_Sigmas_main_Orto('Sigma_i',W_val,1,U_function,V_function,W_Function,z,W_val,U_val,V_val)
    print(abs(Max))
else:
    Draw_3d_Sigmas_main('Sigma_i',W_val,1,U_function,V_function,W_Function,z,W_val,U_val,V_val)


Max_Sigmas_values[0] = Draw_3d_Sigmas('Sigma_x',W_val,1,U_function,V_function,W_Function,z,W_val,U_val,V_val)
Max_Sigmas_values[1] = Draw_3d_Sigmas('Sigma_Y',W_val,2,U_function,V_function,W_Function,z,W_val,U_val,V_val)
Max_Sigmas_values[2] = Draw_3d_Sigmas('Tay_xy',W_val,3,U_function,V_function,W_Function,z,W_val,U_val,V_val)

print("max Sigma_x = ",Max_Sigmas_values[0])
print("max Sigma_Y = ",Max_Sigmas_values[1])
print("max Tay_xy = ",Max_Sigmas_values[2])

Sigma_max_global = ((Max_Sigmas_values[0]**2) + (Max_Sigmas_values[1]**2) - Max_Sigmas_values[0]*Max_Sigmas_values[1] + 3*(Max_Sigmas_values[2]**2))**(1/2)

print("Sigma max global = ",Sigma_max_global)


