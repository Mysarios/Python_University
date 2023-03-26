from sympy import *
import numpy as np
import math as m
import matplotlib as mpl
import matplotlib.pyplot as plt

N=3
L=15
x = Symbol('x')
w_coefs=[]
E=2.1*(10**5)
q_T=1.34/100
h=0.15

#Settings to integral
a = 0
b = L

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
def fuction_q(x,n,l,Coef):
    Pp = Symbol('pi')
    Ll = Symbol('l')
    w = []
    w_2 = []
    w_3 = []
    for i in range(1, n + 1):
        w.append(Coef[i - 1] * sin(m.pi * x * i / Ll))
    return w

#Add w1..n
for i in range(1, N + 1):
    w_coefs.append(Symbol('w'+str(i)))


y = function_w(x,N,L,w_coefs)
q = fuction_q(x,N,L,w_coefs)

y_result=0
q_result=0
for i in range (0,N):
    y_result += y[i]
    q_result += q[i]

y_for_sigma=y_result
y=y_result**2
q=q_result

y_1 = integrate(y,(x,a,b))
q_1 = integrate(q,(x,a,b))

Ee = Symbol('E')
Hh = Symbol('h')
Qq = Symbol('q')
Ll = Symbol('l')
y_1*=Ee*(Hh**3)/24
q_1*=Qq

Es=y_1
dw = []
for i in range (1,N+1):

    buf = y_1.diff(w_coefs[i - 1]) - q_1.diff(w_coefs[i - 1])

    for j in range(1,N+1):
        if j != i :
            buf=buf.subs('w'+str(j),0)

    buf=buf.subs([(Ee,E),(Hh,h),(Qq,q_T),(Ll,L),(pi,m.pi)])
    dw.append(solve(buf)[0])

W_result=0.
# W для конкретного икса
Xx=L/2
pi=m.pi
for i in range(1, N + 1):
    W_result+= dw[i-1] * sin((pi * i * Xx)/L)

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
    for i in range(1, N + 1):
         W_result+= dw[i-1] * sin((pi * i * x_now_1)/L)
    W_graph.append(W_result)
    X_for_graph_1.append(x_now_1)
    x_now_1 +=X_step_1


axs[0].plot(X_for_graph_1,W_graph)

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

axs[1].plot(X_for_graph,Sigma)
plt.show()
#y_for_sigma=y_for_sigma.suba()
