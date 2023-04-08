import sympy as sym
from sympy import *
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

#Выбор типа функции
Choose = 2
Type = 0
#DATA
N=4
E=2.1*(10**5)
q_T=1.34/100
L=12
h=0.12
#graph points
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
                W_result+= dw[(j-1)*New_n + i-1] * (1 - cos(2 * i * B * m.pi / L)) * (1 - cos(2 * j * A * m.pi / L))

    return W_result
#Get Jakobi
def Get_Jacobian(Function, Result_w):
    Jacobian = [0]*N
    Def_Function = [0]*N

    for i in range(0,N):
        Jacobian[i] = [0]*N
        Def_Function[i] = Function.diff('w' + str(i+1))

    for row in range(0,N):
        for column in range(0,N):
            Jacobian[row][column] = Def_Function[row].diff('w' + str(column+1))
            for W_coefs in range(0,N):
                Jacobian[row][column] = Jacobian[row][column].subs('w' + str(W_coefs+1), (Result_w[W_coefs]))

    return Jacobian
#Collect Es
def function_w(x,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w=[]
    w_2=[]
    w_3 = []
    for i in range(1,n+1):
         w.append(Coef[i-1]*sin(m.pi*x*i/Ll))
    for i in range(1, n + 1):
        w_2.append(w[i-1].diff(x))
    for i in range(1, n + 1):
        w_3.append(w_2[i-1].diff(x))
    return w_3
def function_w4(x,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    w_2 = []
    w_3 = []
    for i in range(1, n + 1):
        w.append(Coef[i - 1] * sin(m.pi * x * i / Ll))
    for i in range(1, n + 1):
        w_2.append(w[i - 1].diff(x))
    for i in range(1, n + 1):
        w_3.append(w_2[i - 1].diff(x))

    return w_3
def fuction_q(x,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    for i in range(1, n + 1):
        w.append(Coef[i - 1] * sin(m.pi * x * i / Ll))
    return w
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
                w.append(Coef[(j-1)*New_n + i-1] * (1 - cos(2 * i * y * m.pi / Ll)) * (1 - cos(2 * j * x * m.pi / Ll)))

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
                w.append(Coef[i - 1] * (1 - cos(2 * i * y * m.pi / Ll)) * (1 - cos(2 * j * x * m.pi / Ll)))

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
                w.append(Coef[i - 1] * (1 - cos(2 * i * y * m.pi / Ll)) * (1 - cos(2 * j * x * m.pi / Ll)))

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
        for i in range(1, n + 1):
            w.append(Coef[i - 1] * sin(m.pi * x * i / Ll))
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
                w.append(Coef[i - 1] * (1 - cos(2 * i * y * m.pi / Ll)) * (1 - cos(2 * j * x * m.pi / Ll)))

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
def Lin_Function_Loop(w_coefs,y_1,q_1,fig, axs):
    Ee = Symbol('E')
    Hh = Symbol('h')
    Qq = Symbol('q')
    Ll = Symbol('l')
    W_Result = [0] * (Size + 2)

    W_graph = []
    Q_now = 0.00001
    Q_step = q_T/(Size/Max_q_T)
    q_for_graph = []
    # Цикл по разным q
    for j in range(1, Size + 2):
        dw = []
        print("q_now_Linear = ", Q_now)
        for i in range(1, N + 1):
            buf = y_1.diff(w_coefs[i - 1]) - q_1.diff(w_coefs[i - 1])
            # Сокращение ~~0 коэф.
            for j in range(1, N + 1):
                if j != i:
                    buf = buf.subs('w' + str(j), 0)
            buf = buf.subs([(Ee, E), (Hh, h), (Qq, Q_now), (Ll, L), (pi, m.pi)])
            dw.append(solve(buf)[0])
        W_graph.append(Get_W(L / 2,dw))
        print(Get_W(L / 2,dw))
        Q_now += Q_step
        q_for_graph.append(Q_now)
        W_Result.append(dw)

    axs.plot(W_graph, q_for_graph)
    plt.show()

    return W_Result
def Lin_Function(w_coefs,y_1,q_1,Q__T):
    Ee = Symbol('E')
    Hh = Symbol('h')
    Qq = Symbol('q')
    Ll = Symbol('l')
    dw = []

    for i in range(1, N + 1):
        buf = y_1.diff(w_coefs[i - 1]) - q_1.diff(w_coefs[i - 1])
        # Сокращение ~~0 коэф.
        for j in range(1, N + 1):
              if j != i:
                  buf = buf.subs('w' + str(j), 0)
        buf = buf.subs([(Ee, E), (Hh, h), (Qq, Q__T), (Ll, L), (pi, m.pi)])
        dw.append(solve(buf)[0])

    return dw
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
#Не линейный варик
def Ne_Lin_Function_Loop(w_coefs,y_1,q_1,yn_1,fig, axs):
    Ee = Symbol('E')
    Hh = Symbol('h')
    Qq = Symbol('q')
    Ll = Symbol('l')
    W_Result = [0]*(Size+2)
    W_Result[0] = [0] * N
    Q_now = 0.00001
    Q_step = q_T/(Size/Max_q_T)
    q_for_graph =[]
    #Цикл по разным q

    q_for_graph.append(Q_now)
    for j in range(1,Size+2):
        W_Result[j]=W_Result[j-1]
        print("q_now_Not_Linear = ", Q_now)


        Buf_Function = y_1 + yn_1 - q_1
        Main_Function = Buf_Function.subs([(Ee, E), (Hh, h), (Qq, Q_now), (Ll, L), (pi, m.pi)])
        W_Result[j] = Nuton_Iter(Main_Function,eps, W_Result[j], w_coefs)

        Q_now += Q_step
        q_for_graph.append(Q_now)
        print(Get_W(L / 2,W_Result[j]))

    W_graph = []
    for i in range(0,Size+2):
        W_graph.append(Get_W(L/2,W_Result[i]))

    axs.plot(W_graph,q_for_graph)
    return W_Result
def Ne_Lin_Function(w_coefs,y_1,q_1,yn_1,Q__T):
    Ee = Symbol('E')
    Hh = Symbol('h')
    Qq = Symbol('q')
    Ll = Symbol('l')

    W_Result = [0]*N

    Buf_Function = y_1 + yn_1 - q_1
    Main_Function = Buf_Function.subs([(Ee, E), (Hh, h), (Qq, Q__T), (Ll, L), (pi, m.pi)])

    W_Result = Nuton_Iter(Main_Function, eps, W_Result, w_coefs)

    return W_Result
#Nut iter
def Nuton_Iter(F,eps,w0,w_coefs):
    Now_eps = 1
    All_Results = []
    Res_now = []
    Res_now = w0
    Res_Last_now = Res_now
    Res_new = []
    Count_Iter = 0
    Check_Loop =1
    Loop = 0

    New_Eps = [0]*N

    while(Now_eps > eps):
        Count_Iter += 1
        Max_eps = 0

        Jacobi = Get_Jacobian(F,Res_now)
        Jackobi_Matrix = sym.Matrix(Jacobi)
        Jacobi_Invariant = Jackobi_Matrix.inv()

        Res_new = Get_New_iterarion(F,Jacobi_Invariant,Res_now,Check_Loop)
        np.array(Res_new).astype(np.float64)

        for i in range(0,N):
            New_Eps[i] = abs(Res_new[i] - Res_now[i])
            if(New_Eps[i] > Max_eps):
                Max_eps = New_Eps[i]

        Res_Last_now = Res_now
        Now_eps = Max_eps
        Res_now = Res_new
        All_Results.append(Res_now)

        if Count_Iter > 10:
            Check_Loop = Check_Loop/10
            Res_now = Res_Last_now
            if All_Results[Count_Iter-3][0] < All_Results[Count_Iter-2][0] < All_Results[Count_Iter-1][0]:
                if Check_Loop != 1:
                    Check_Loop *= 10
                Res_now = Res_new
            if All_Results[Count_Iter-3][0] > All_Results[Count_Iter-2][0] > All_Results[Count_Iter-1][0]:
                if Check_Loop != 1:
                    Check_Loop *= 10
                Res_now = Res_new

        if Count_Iter > 50:
            Now_eps = 0

    return Res_now

#Создание апрк. функции
if Choose == 1 or Choose == 0:
    w_coefs = Create_w()
else:
    w_coefs = Create_w_Plane()
#print(w_coefs)

#Создание частей Es
if Choose == 1 or Choose == 0:
    y = function_w(x,N,L,w_coefs)
    q = fuction_q(x,N,L,w_coefs)
    y_neLin = function_w4(x,N,L,w_coefs)
else:
    Es_1 = Es_1_function(x,Yy,N,L,w_coefs)
    Es_2 = Es_2_function(x,Yy,N,L,w_coefs)
    Es_3 = Es_3_function(x,Yy,N,L,w_coefs)
    Es_4 = Es_4_function(x,Yy,N,L,w_coefs)
    Es_5 = Es_5_function(x,Yy,N,L,w_coefs)


#Соед.Частей Es
if Choose == 1 or Choose == 0:
    y_result=0
    q_result=0
    y_neLin_result=0
    for i in range (0,N):
        y_result += y[i]
        q_result += q[i]
        y_neLin_result += y_neLin[i]
    y_for_sigma=y_result
else:
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
if Choose == 1 or Choose == 0:
    y=y_result**2
    y_neLin= y_neLin_result**4
    q=q_result
    y_1 = integrate(y,(x,a,b))
    q_1 = integrate(q,(x,a,b))
    y_nelin_1 = integrate(y_neLin,(x,a,b))
else:
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

if Choose == 1 or Choose == 0:
    y_1*=Ee*(Hh**3)/24
    q_1*=Qq
    y_nelin_1*=(-1)*(Ee*(Hh**5)*Mm/120)
else:
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
fig,axs = plt.subplots(2)
x_now_1=0
X_count_1 = L/0.5
X_step_1 = L/X_count_1
X_for_graph_1 =[]
W_graph =[]

print("Hi")
if(Choose == 0):
    dw = Lin_Function(w_coefs,y_1,q_1,Q__T)
    dw_2 = Ne_Lin_Function(w_coefs, y_1, q_1, y_nelin_1,Q__T)

    print(Get_W(L / 2, dw))
    print(Get_W(L / 2, dw_2))

    # W для всех иксов
    if Choose == 0:
        for i in range(1, int(X_count_1) + 2):
            W_result = 0
            W_result = Get_W(x_now_1, dw)
            W_graph.append(W_result)
            X_for_graph_1.append(x_now_1)
            x_now_1 += X_step_1
        # Вывод графика

        axs[0].plot(X_for_graph_1, W_graph)

        # Просчет сигм
        y_for_sigma = y_for_sigma.subs([(Ll, L), (pi, m.pi)])
        Sigma = []
        z = h / 2
        x_now = 0
        X_count = L / 0.5
        X_step = L / X_count
        X_for_graph = []
        for i in range(1, int(X_count) + 2):
            X_for_graph.append(x_now)
            E_x = y_for_sigma.subs(x, x_now)
            x_now += X_step
            for j in range(1, N + 1):
                E_x = E_x.subs('w' + str(j), dw[j - 1])

            buf_sigma = -1 * E * z * E_x
            Sigma.append(buf_sigma)

        # Вывод сигм
        axs[1].plot(X_for_graph, Sigma)
        plt.show()
if (Choose == 1):
    fig, axs = plt.subplots(1)
    dw= Ne_Lin_Function_Loop(w_coefs,y_1,q_1,y_nelin_1,fig, axs)
    dw = Lin_Function_Loop(w_coefs, y_1, q_1,fig, axs)
if (Choose == 2):
    dw = Lin_Function_Plane(w_coefs,Es_1,Es_2,Es_3,Es_4,Es_5,N)
    Result = Get_W_Plane(A/2,B/2, dw)
    print(Result)
    print(dw)



