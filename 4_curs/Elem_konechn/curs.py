from sympy import *
import numpy as np
import math as m
import matplotlib.pyplot as plt

Close_count = 3
Size_graph = 100

Ww = Symbol('w')
Qq = Symbol('q')
Ee = Symbol('E')
Ll = Symbol('l')
Ii = Symbol('I')
Xx = Symbol('x')
PI = Symbol('pi')

I = 317.75
E = 120 * 1000
q = 10
l = 5
def Get_analyt_result():
    print("\nStart Analyt\n")
    C_coefs = [Symbol('c1'),Symbol('c2'),Symbol('c3'),Symbol('c4')]
    Get_Esp = 0
    Get_Esp += (1/24)*q/(E*I)*(Xx**4)
    Get_Esp += (1/6)*(Xx**3)*C_coefs[0]
    Get_Esp += (1/2)*(Xx**2)*C_coefs[1]
    Get_Esp += Xx*C_coefs[2]
    Get_Esp += C_coefs[3]

    print("Analyt W(x) =",Get_Esp)

    Dw = [Get_Esp]
    for i in range(0, 4):
        Dw.append(Dw[i].diff(Xx))

    print("Dw[1]", Dw[1])
    print("Dw[2]", Dw[2])
    print("Dw[3]", Dw[3])
    print("Dw[4]", Dw[4])

    L = E * I * Dw[4] - q
    print("L = ",L)
    System_get = [0] * 4
    System_get[0] = Dw[0].subs([(Xx, 0),])
    System_get[1] = Dw[0].subs([(Xx, l),])
    System_get[2] = Dw[1].subs([(Xx, 0),])
    System_get[3] = Dw[1].subs([(Xx, l),])

    for i in range(0, 4):
        print("System_get[", i + 1, "] =", System_get[i])

    print("C_coefs =", C_coefs)
    Result = 0
    for solution in linsolve(System_get, C_coefs):
        Result = solution

    print("Result_coefs =",Result)

    Result_y = []

    Result_num = []
    for res in Result:
        Result_num.append(res.subs([(Ll, l), (Ii, I)]))

    print("Coefs num =", Result_num)

    Result_W_x = Get_Esp
    Result_W_x = Result_W_x.subs([(Ll, l)])
    for i in range(0, len(C_coefs)):
        Result_W_x = Result_W_x.subs([(C_coefs[i], Result_num[i])])

    print("W(x) = ", Result_W_x)

    x = np.linspace(0, l, Size_graph)
    y_W = [0] * Size_graph

    for i in range(0, Size_graph):
        y_W[i] = Result_W_x.subs([(Xx, x[i])])

    Result_y.append(y_W)
    # plt.plot(x, y_W)
    # plt.show()

    Result_M_x = -Ee * Ii * Dw[2]
    Result_M_x = Result_M_x.subs([(Ll, l), (Ee, E), (Ii, I)])
    for i in range(0, len(C_coefs)):
        Result_M_x = Result_M_x.subs([(C_coefs[i], Result_num[i])])
    print("M(x) = ", Result_M_x)

    y_M = [0] * Size_graph

    for i in range(0, Size_graph):
        y_M[i] = Result_M_x.subs([(Xx, x[i])])
    Result_y.append(y_M)


    Result_Q_x = Result_M_x.diff(Xx)
    print("Q(x) = ", Result_Q_x)

    y_Q = [0] * Size_graph

    for i in range(0, Size_graph):
        y_Q[i] = Result_Q_x.subs([(Xx, x[i])])

    Result_y.append(y_Q)
    # plt.plot(x, y_Q)
    # plt.show()

    return Result_y
def Get_Bubnov_result(Count_iter):
    print("\nStart Bubnov(",Count_iter,")\n")

    Fi_x = []
    for i in range(0, Count_iter):
        Fi_x.append(GetSin(i + 1))

    W_coefs = Get_w_coefs(Count_iter)
    Get_Esp = 0
    for i in range(0, Count_iter):
        Get_Esp += get_W_sin(W_coefs[i], i + 1)

    print("Analyt W(x) =", Get_Esp)
    print("Fi =", Fi_x)

    Dw = [Get_Esp]
    for i in range(0, 4):
        Dw.append(Dw[i].diff(Xx))

    print("Dw[1]", Dw[1])
    print("Dw[2]", Dw[2])
    print("Dw[3]", Dw[3])
    print("Dw[4]", Dw[4])

    L = Ee * Ii * Dw[4] - Qq

    Integ = [0] * Count_iter
    for i in range(0, Count_iter):
        Integ[i] = Fi_x[i] * L

    for i in range(0, Count_iter):
        print("Integ[", i + 1, "] =", Integ[i])

    Get_Integ = [0] * Count_iter
    for i in range(0, Count_iter):
        Get_Integ[i] = integrate(Integ[i], (Xx, 0, Ll))

    for i in range(0, Count_iter):
        print("Get_Integ[", i + 1, "] =", Get_Integ[i])

    for i in range(0, len(Get_Integ)):
        print("I[", i, "] =", Get_Integ[i])
    print("All I =", Get_Integ)
    print("W_c =", W_coefs)
    for solution in linsolve(Get_Integ, W_coefs):
        Result = solution

    Result_y = []
    Result_num = []
    for res in Result:
        Result_num.append(res.subs([(Ee, E), (Qq, q), (Ll, l), (Ii, I)]))

    print("Coefs num =", Result_num)

    Result_W_x = Get_Esp
    Result_W_x = Result_W_x.subs([(Ll, l)])
    for i in range(0, Count_iter):
        Result_W_x = Result_W_x.subs([(W_coefs[i], Result_num[i])])

    print("W(x) = ", Result_W_x)

    x = np.linspace(0, l, Size_graph)
    y_W = [0] * Size_graph

    for i in range(0, Size_graph):
        y_W[i] = Result_W_x.subs([(Xx, x[i])])

    Result_y.append(y_W)
    Result_M_x = -Ee * Ii * Dw[2]
    Result_M_x = Result_M_x.subs([(Ll, l), (Ee, E), (Ii, I)])
    for i in range(0, Count_iter):
        Result_M_x = Result_M_x.subs([(W_coefs[i], Result_num[i])])

    print("M(x) = ", Result_M_x)

    y_M = [0] * Size_graph

    for i in range(0, Size_graph):
        y_M[i] = Result_M_x.subs([(Xx, x[i])])

    Result_y.append(y_M)
    Result_Q_x = Result_M_x.diff(Xx)

    print("Q(x) = ", Result_Q_x)

    y_Q = [0] * Size_graph

    for i in range(0, Size_graph):
        y_Q[i] = Result_Q_x.subs([(Xx, x[i])])

    Result_y.append(y_Q)
    return Result_y

def Get_w_coefs(Count):
    w_coefs = []
    for i in range(1, Count + 1):
        w_coefs.append(Symbol('w' + str(i)))
    return w_coefs
def GetSin(j):
    return (sin((2*j-1)*m.pi*Xx/ Ll))**2
def get_W_sin(Wc,j):
    return Wc*(sin((2 * j - 1) * m.pi * Xx / Ll)) ** 2

x_graph = np.linspace(0, l, Size_graph)
y_array = [0]*4
for i in range(0,3):
    y_array[i] = (Get_Bubnov_result(i+1))
y_array[3] = (Get_analyt_result())

for i in range(0,4):
    plt.plot(x_graph, y_array[i][0])
plt.legend([f"W(x) при N = 1" , f"W(x) при N = 2", f"W(x) при N = 3", f"W(x) analyt"])
plt.show()

for i in range(0,4):
    plt.plot(x_graph, y_array[i][1])
plt.legend([f"M(x) при N = 1" , f"M(x) при N = 2", f"M(x) при N = 3", f"M(x) analyt"])
plt.show()

for i in range(0,4):
    plt.plot(x_graph, y_array[i][2])
plt.legend([f"Q(x) при N = 1" , f"Q(x) при N = 2", f"Q(x) при N = 3", f"Q(x) analyt"])
plt.show()



