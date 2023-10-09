from sympy import *
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

N=1
L=15
x = Symbol('x')
E=2.1*(10**5)
q_T=1.34/100
h=0.15

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
#return W
def Get_W(X,dw):
    W_result = 0.
    Xx=L/2
    pi=m.pi
    for i in range(1, N + 1):
        W_result+= dw[i-1] * sin((pi * i * Xx)/L)
    return W_result
#Diff w to Es
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
    w_4 = []
    w_5 = []
    for i in range(1, n + 1):
        w.append(Coef[i - 1] * sin(m.pi * x * i / Ll))
    for i in range(1, n + 1):
        w_2.append(w[i - 1].diff(x))
    for i in range(1, n + 1):
        w_3.append(w_2[i - 1].diff(x))
    for i in range(1, n + 1):
        w_4.append(w_3[i - 1].diff(x))
    for i in range(1, n + 1):
        w_5.append(w_4[i - 1].diff(x))
    return w_5
#Collect q to Es
def fuction_q(x,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    w_2 = []
    w_3 = []
    for i in range(1, n + 1):
        w.append(Coef[i - 1] * sin(m.pi * x * i / Ll))
    return w
#Loop for Nuton

#Линейный варик
def Lin_Function(w_coefs):
    W_Result = 0







    return W_Result
#Не линейный варик
def Ne_Lin_Function(w_coefs):
    W_Result=0
    return W_Result
#Nut iter
def Nuton_Iter(F,dF,eps,w0):
    Now_eps = 0
    Res_w = [w0,0]
    F_buf = 0
    Df_buf = 0

    while(Now_eps > eps):
        E_x_1 = 0
        E_x_2 = 0
        for j in range(1, N + 1):
            E_x_1 = F.subs('w' + str(j), dw[j - 1])
            E_x_2 = dF.subs('w' + str(j), dw[j - 1])

        F_buf = E_x_1
        Df_buf = E_x_2
        Now_eps = abs(F_buf/Df_buf)
        Res_w[0]-=F_buf/Df_buf
    for j in range(1, N + 1):
        E_x_1 = F.subs('w' + str(j), dw[j - 1])
    Res_w[1] = E_x_1

    return Res_w

Choose =0
w_coefs = Create_w()

y = function_w(x,N,L,w_coefs)
q = fuction_q(x,N,L,w_coefs)
if(Choose == 1)
    y_neLin = function_w4(x,N,L,w_coefs)

y_result=0
q_result=0
y_neLin_result=0
for i in range (0,N):
    y_result += y[i]
    q_result += q[i]
    if (Choose == 1)
        y_neLin_result += y_neLin[i]

y_for_sigma=y_result
y=y_result**2
y_neLin= y_neLin_result**4
q=q_result

y_1 = integrate(y,(x,a,b))
q_1 = integrate(q,(x,a,b))
y_nelin_1 = integrate(y_neLin,(x,a,b))

Ee = Symbol('E')
Hh = Symbol('h')
Qq = Symbol('q')
Ll = Symbol('l')
y_1*=Ee*(Hh**3)/24
q_1*=Qq
y_nelin_1*=Ee*(Hh**5)*m/120

Es=y_1
dw = []
for i in range (1,N+1):
    buf = y_1.diff(w_coefs[i - 1]) - q_1.diff(w_coefs[i - 1])
    for j in range(1,N+1):
        if j != i :
            buf=buf.subs('w'+str(j),0)

    buf=buf.subs([(Ee,E),(Hh,h),(Qq,q_T),(Ll,L),(pi,m.pi)])
    dw.append(solve(buf)[0])

# W для конкретного икса
W_result=0.
X_for_W=L/2
W_result = Get_W(X_for_W,dw)
print(W_result)

# W для всех иксов
fig,axs = plt.subplots(2)
x_now_1=0
X_count_1 = L/0.5
X_step_1 = L/X_count_1
X_for_graph_1 =[]
W_graph =[]
for i in range (1,int(X_count_1)+2):
    W_result=0
    W_result = Get_W(x_now_1,dw)
    W_graph.append(W_result)
    X_for_graph_1.append(x_now_1)
    x_now_1 +=X_step_1
#Вывод графика
axs[0].plot(X_for_graph_1,W_graph)


#Просчет сигм
y_for_sigma=y_for_sigma.subs([(Ll,L),(pi,m.pi)])
Sigma =[]
z=h/2
x_now=0
X_count = L/0.5
X_step = L/X_count
X_for_graph=[]
for i in range (1,int(X_count)+2):
    X_for_graph.append(x_now)
    E_x =y_for_sigma.subs(x,x_now)
    x_now+=X_step
    for j in range(1,N+1):
            E_x=E_x.subs('w'+str(j),dw[j-1])

    buf_sigma = -1*E*z*E_x
    Sigma.append(buf_sigma)

#Вывод сигм
axs[1].plot(X_for_graph,Sigma)
plt.show()



