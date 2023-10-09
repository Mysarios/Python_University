import sympy as sym
from sympy import *
import numpy as np
import math as m
import matplotlib as mpl
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

#Выбор типа функции
Choose = 2
Type = 2
#DATA
N=4
E=2.1*(10**5)
q_T=1.34/100
L=15
h=0.15
z_val = h/2
#graph points
size_Graph = 20
Size = 30
Max_q_T = 2
A=L
B=L
x = Symbol('x')
Yy = Symbol('y')
Q__T = q_T
Mm=10**5
eps = 0.001
nu = 0.3
#Settings to integral
a = 0
b = L

#Создать w функцию
def Draw_3d_W(dw):
    x_array = []
    y_array = []
    z_array = [0]*Size
    for i in range (0,Size):
        z_array[i] = [0]*Size
    step = A/(Size-1)
    for i in range(0,Size):
        x_array.append(i*step)
        y_array.append(i*step)

    for i in range(0,Size):
        for j in range(0,Size):
            z_array[i][j] = (Get_W_Plane(x_array[i],y_array[j],dw))

    z_array=np.array(z_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, z_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('W (x,y)')
    ax.set_title('W (x,y)')
    plt.show()
def Sigmas(dw,Coefs):
    w =Get_W_To_Sigmas(Coefs)
    x = Symbol('x')
    y = Symbol('y')
    All_w = 0
    for i in range(0,N):
        All_w += w[i]

    w_for_Sigma_x = All_w.diff(x)
    w_for_Sigma_x = w_for_Sigma_x.diff(x)

    w_for_Sigma_y = All_w.diff(y)
    w_for_Sigma_y = w_for_Sigma_y.diff(y)

    w_for_tay = All_w.diff(x)
    w_for_tay = w_for_tay.diff(y)

    x_array = []
    y_array = []
    Sigma_X = [0] * Size
    Sigma_Y = [0] * Size
    Tay = [0] * Size

    for i in range(0, Size):
        Sigma_X[i] = [0] * Size
        Sigma_Y[i] = [0] * Size
        Tay[i] = [0] * Size

    step = A / (Size - 1)
    for i in range(0, Size):
        x_array.append(i * step)
        y_array.append(i * step)


    Plot_Sigmas_X(Sigma_X,x_array,y_array,w_for_Sigma_x,dw,w_for_Sigma_x,w_for_Sigma_y,w_for_tay)
    Plot_Sigmas_Y(Sigma_X, x_array, y_array, w_for_Sigma_y, dw,w_for_Sigma_x,w_for_Sigma_y,w_for_tay)
    Plot_Tay(Sigma_X, x_array, y_array, w_for_tay, dw,w_for_Sigma_x,w_for_Sigma_y,w_for_tay)

def Plot_Sigmas_X(Sigmas_array,x_array,y_array,dw_Sigmas,dw,dw_X,dw_Y,dw_xy):
    y = Symbol('y')
    x = Symbol('x')
    max = 0;
    for W_coefs in range(0, N):
        #dw_Sigmas = dw_Sigmas.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_X = dw_X.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_Y = dw_Y.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_xy = dw_xy.subs('w' + str(W_coefs + 1), (dw[W_coefs]))

    for i in range(0,Size):
        for j in range(0,Size):
            Sigmas_array[i][j] = (-z_val*E/(1 - nu**2)*(dw_X + nu*dw_Y)).subs([(x,x_array[i]),(y,y_array[j])])
            if Sigmas_array[i][j] > max:
                max =Sigmas_array[i][j]

    print("Max Sigma x =",max)

    Sigmas_array = np.array(Sigmas_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Sigmas_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Sigma (x,y)')
    ax.set_title('Sigma X')
    #ax.xlabel('x')
    #ax.ylabel('y')
    #ax.zlabel('Sigma(x,y)')
    plt.show()
def Plot_Sigmas_Y(Sigmas_array, x_array, y_array, dw_Sigmas, dw,dw_X,dw_Y,dw_xy):
    y = Symbol('y')
    x = Symbol('x')
    max = 0;
    for W_coefs in range(0, N):
        # dw_Sigmas = dw_Sigmas.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_X = dw_X.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_Y = dw_Y.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_xy = dw_xy.subs('w' + str(W_coefs + 1), (dw[W_coefs]))

    for i in range(0, Size):
        for j in range(0, Size):
            Sigmas_array[i][j] = (-z_val*E/(1 - nu**2)*(dw_Y + nu*dw_X)).subs([(x, x_array[i]), (y, y_array[j])])
            if Sigmas_array[i][j] > max:
                max =Sigmas_array[i][j]

    print("Max Sigma y =", max)
    Sigmas_array = np.array(Sigmas_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Sigmas_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Sigma (x,y)')
    ax.set_title('Sigma Y')
    plt.show()
def Plot_Tay(Sigmas_array, x_array, y_array, dw_Sigmas, dw,dw_X,dw_Y,dw_xy):
    y = Symbol('y')
    x = Symbol('x')
    max = 0;

    for W_coefs in range(0, N):
        # dw_Sigmas = dw_Sigmas.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_X = dw_X.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_Y = dw_Y.subs('w' + str(W_coefs + 1), (dw[W_coefs]))
        dw_xy = dw_xy.subs('w' + str(W_coefs + 1), (dw[W_coefs]))

    for i in range(0, Size):
        for j in range(0, Size):
            Sigmas_array[i][j] = (-z_val*E/(1 + nu)*dw_xy).subs([(x, x_array[i]), (y, y_array[j])])
            if Sigmas_array[i][j] > max:
                max =Sigmas_array[i][j]

    print("Max tay =", max)
    Sigmas_array = np.array(Sigmas_array)
    x_array = np.array(x_array)
    y_array = np.array(y_array)
    X, Y = np.meshgrid(x_array, y_array)
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.plot_surface(X, Y, Sigmas_array, cmap='viridis', edgecolor='green')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('Tay (x,y)')
    ax.set_title('Tay')
    plt.show()
def Create_w():
    w_coefs = []
    # Add w1..n
    for i in range(1, N + 1):
        w_coefs.append(Symbol('w' + str(i)))
    return w_coefs
def Create_w_Plane():
    w_coefs = []
    # Add w1..n
    #Hard
    for i in range(1, N + 1):
        w_coefs.append(Symbol('w' + str(i)))
    return w_coefs

#return W
def Get_W(X,dw):
    W_result = 0.
    for i in range(1, N + 1):
        W_result+= dw[i-1] * sin((m.pi * i * X)/L)

    return W_result
def Get_W_To_Sigmas(Coef):
    w = []
    y=Symbol('y')
    x=Symbol('x')
    New_n = N ** (1 / 2.0)
    New_n = int(New_n)
    # Hard
    if Type == 0:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j - 1) * New_n + i - 1] * (1 - cos(2 * i * y * m.pi / L)) * (1 - cos(2 * j * x * m.pi / L)))
    if Type == 1:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j - 1) * New_n + i - 1] * sin(i * y * m.pi / L) * sin(j * x * m.pi / L))
    if Type == 2:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j - 1) * New_n + i - 1] * (1 - cos(2 * i * y * m.pi / L)) * sin(j * x * m.pi / L))

    return w
def Get_W_Plane(X_val,Y_val,dw):
    W_result = 0.
    New_n = N ** (1 / 2.0)
    New_n = int(New_n)


    # Hard
    if Type == 0:
        for j in range(1, New_n+1):
            for i in range(1, New_n+1):
                W_result+= dw[(j-1)*New_n + i-1] * (1 - cos(2 * i * Y_val * m.pi / A)) * (1 - cos(2 * j * X_val * m.pi / B))

    if Type == 1:
        for j in range(1, New_n+1):
            for i in range(1, New_n+1):
                W_result+= dw[(j-1)*New_n + i-1] * sin(i * Y_val * m.pi / L) * sin(j * X_val * m.pi / L)

    if Type == 2:
        for j in range(1, New_n+1):
            for i in range(1, New_n+1):
                W_result+= dw[(j-1)*New_n + i-1] * (1 - cos(2 * i * Y_val * m.pi / L)) * sin(j * X_val * m.pi / L)

    return W_result

def Es_1_function(x,y,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    w_2 = []
    w_3 = []
    New_n = N ** (1 / 2.0)
    New_n = int(New_n)

    # Hard
    if Type == 0:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * (1 - cos(2 * i * y * m.pi / L)) * (1 - cos(2 * j * x * m.pi / L)))

    if Type == 1:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * sin(i * y * m.pi / L) * sin(j * x * m.pi / L))
    if Type == 2:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * (1 - cos(2 * i * y * m.pi / Ll)) * sin(j * x * m.pi / L))

    for i in range(1, n + 1):
        w_2.append(w[i - 1].diff(x))
    for i in range(1, n + 1):
        w_3.append(w_2[i - 1].diff(x))

    return w_3
def Es_2_function(x,y,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    w_2_x = []
    w_2_y = []
    w_3_y = []
    w_3_x = []
    w_res = []
    New_n = N ** (1 / 2.0)
    New_n = int(New_n)
    # Hard
    if Type == 0:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] *(1 - cos(2 * i * y * m.pi / L)) * (1 - cos(2 * j * x * m.pi / L)))
    if Type == 1:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * sin(i * y * m.pi / L) * sin(j * x * m.pi / L))
    if Type == 2:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[i - 1] * (1 - cos(2 * i * y * m.pi / Ll)) * sin(j * x * m.pi / L))

    for i in range(1, n + 1):
        w_2_x.append(w[i - 1].diff(x))
        w_2_y.append(w[i - 1].diff(y))
    for i in range(1, n + 1):
        w_3_y.append(w_2_y[i - 1].diff(y))
        w_3_x.append(w_2_x[i - 1].diff(x))
    for i in range(1, n + 1):
        w_res.append(w_3_y[i-1] * w_3_x[i-1])
    return w_res
def Es_3_function(x,y,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    w_2 = []
    w_3 = []
    New_n = N**(1/2.0)
    New_n = int(New_n)
    #Hard
    if Type == 0:
        for j in range(1,New_n+1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * ( 1 - cos(2*i*y*m.pi/L)) * ( 1 - cos(2*j*x*m.pi/L)) )
    if Type == 1:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * sin(i * y * m.pi / L) * sin(j * x * m.pi / L))
    if Type == 2:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[i - 1] * (1 - cos(2 * i * y * m.pi / Ll)) * sin(j * x * m.pi / L))

    for i in range(1, n + 1):
        w_2.append(w[i - 1].diff(y))
    for i in range(1, n + 1):
        w_3.append(w_2[i - 1].diff(y))
    return w_3
def Es_4_function(x,y,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    w_2 = []
    w_3 = []
    New_n = N ** (1 / 2.0)
    New_n = int(New_n)
    # Hard
    if Type == 0:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * (1 - cos(2 * i * y * m.pi / L)) * (1 - cos(2 * j * x * m.pi / L)))
    if Type == 1:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * sin(i * y * m.pi / L) * sin(j * x * m.pi / L))
    if Type == 2:
        for j in range(1, New_n + 1):
            for i in range(1, n + 1):
                w.append(Coef[i - 1] * (1 - cos(2 * i * y * m.pi / L)) * sin(j * x * m.pi / L))
    for i in range(1, n + 1):
        w_2.append(w[i - 1].diff(x))
    for i in range(1, n + 1):
        w_3.append(w_2[i - 1].diff(y))
    return w_3
def Es_5_function(x, y, n, l, Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    New_n = n ** (1 / 2.0)
    New_n = int(New_n)
    # Hard
    if Type == 0:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * (1 - cos(2 * i * y * m.pi / L)) * (1 - cos(2 * j * x * m.pi / L)))

    if Type == 1:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[(j-1)*New_n + i-1] * sin(i * y * m.pi / Ll) * sin(j * x * m.pi / Ll))

    if Type == 2:
        for j in range(1, New_n + 1):
            for i in range(1, New_n + 1):
                w.append(Coef[i - 1] * (1 - cos(2 * i * y * m.pi / Ll)) * sin(j * x * m.pi / Ll))

    return w
#Loop for Nuton
def Get_New_iterarion(Function,Jackobi_inv,W,a):
    Def_Function = [0] * N

    for i in range(0, N):
        Def_Function[i] = Function.diff('w' + str(i + 1))
        for W_coefs in range(0, N):
            Def_Function[i] = Def_Function[i].subs('w' + str(W_coefs + 1), (W[W_coefs]))

    Def_Function = sym.Matrix(Def_Function)
    W = sym.Matrix(W)
    W = W - a*(Jackobi_inv * Def_Function)
    return W
#Линейный варик
def Lin_Function_Plane(w_coefs,Es_1,Es_2,Es_3,Es_4,Es_5,N):
    Ee = Symbol('E')
    Hh = Symbol('h')
    Qq = Symbol('q')
    Ll = Symbol('l')
    dw = []

    for i in range(1, N + 1):
        buf = Es_1.diff(w_coefs[i - 1]) + Es_2.diff(w_coefs[i - 1]) + Es_3.diff(w_coefs[i - 1]) + Es_4.diff(w_coefs[i - 1]) - Es_5.diff(w_coefs[i - 1])
        #Сокращение ~~0 коэф.
        for j in range(1, N + 1):
              if j != i:
                  buf = buf.subs('w' + str(j), 0)
        buf = buf.subs([(Ee, E), (Hh, h), (Qq, Q__T), (Ll, L), (pi, m.pi)])
        dw.append(solve(buf)[0])

    return dw

#Создание апрк. функции
w_coefs = Create_w_Plane()

#Создание частей Es
Es_1 = Es_1_function(x,Yy,N,L,w_coefs)
Es_2 = Es_2_function(x,Yy,N,L,w_coefs)
Es_3 = Es_3_function(x,Yy,N,L,w_coefs)
Es_4 = Es_4_function(x,Yy,N,L,w_coefs)
Es_5 = Es_5_function(x,Yy,N,L,w_coefs)


#Соед.Частей Es

Es_1_result = 0
Es_2_result = 0
Es_3_result = 0
Es_4_result = 0
Es_5_result = 0
for i in range (0,N):
    Es_1_result += Es_1[i]
    Es_2_result += Es_2[i]
    Es_3_result += Es_3[i]
    Es_4_result += Es_4[i]
    Es_5_result += Es_5[i]


#Приведение к нужному виду
Es_1 = Es_1_result**2
Es_2 = Es_2_result
Es_3 = Es_3_result**2
Es_4 = Es_4_result**2
Es_5 = Es_5_result

Es_1 = integrate(Es_1, (x, a, A))
Es_2 = integrate(Es_2, (x, a, A))
Es_3 = integrate(Es_3, (x, a, A))
Es_4 = integrate(Es_4, (x, a, A))
Es_5 = integrate(Es_5, (x, a, A))

Es_1 = integrate(Es_1, (Yy, a, B))
Es_2 = integrate(Es_2, (Yy, a, B))
Es_3 = integrate(Es_3, (Yy, a, B))
Es_4 = integrate(Es_4, (Yy, a, B))
Es_5 = integrate(Es_5, (Yy, a, B))

#Создание символов/Коэф.Интегралов
Ee = Symbol('E')
Hh = Symbol('h')
Qq = Symbol('q')
Ll = Symbol('l')
D = ( Ee*(Hh**3) )/ (12 * (1 - nu**2) )

Es_1 *= D/2
Es_2 *= nu*D
Es_3 *= D/2
Es_4 *= (1 - nu)*D
Es_5 *= Qq

#Массив производных
#Es=y_1
dw = []
dw_2 = []

#Данные для графика
x_now_1=0
X_count_1 = L/0.5
X_step_1 = L/X_count_1
X_for_graph_1 =[]
W_graph =[]

print("Hi")


dw = Lin_Function_Plane(w_coefs,Es_1,Es_2,Es_3,Es_4,Es_5,N)
Result = Get_W_Plane(A/2,B/2, dw)
print(Result)
print(dw)
print(w_coefs)
Draw_3d_W(dw)
Sigmas(dw,w_coefs)




