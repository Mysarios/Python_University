from sympy import *
import numpy as np
import math as m
import matplotlib.pyplot as plt

Close_count = 2
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

def Get_w_coefs():
    w_coefs = []
    for i in range(1, Close_count + 1):
        w_coefs.append(Symbol('w' + str(i)))
    return w_coefs
def GetSin(j):
    return (sin((2*j-1)*m.pi*Xx/ Ll))**2
def get_W_sin(Wc,j):
    return Wc*(sin((2 * j - 1) * m.pi * Xx / Ll)) ** 2

Fi_x = []
for i in range(0,Close_count):
    Fi_x.append(GetSin(i+1))

W_coefs = Get_w_coefs()

Get_Esp = 0
for i in range(0,Close_count):
    Get_Esp += get_W_sin(W_coefs[i],i+1)

print("Esp=",Get_Esp)
print("Fi =",Fi_x)

Dw = [Get_Esp]

for i in range(0,4):
    Dw.append(Dw[i].diff(Xx))

print("Dw[4]",Dw[1])
print("Dw[4]",Dw[2])
print("Dw[4]",Dw[3])
print("Dw[4]",Dw[4])

L = Ee*Ii*Dw[4]-Qq

Integ = [0]*Close_count
for i in range(0,Close_count):
    Integ[i] = Fi_x[i]*L
for i in range(0,Close_count):
    print("Integ[",i+1,"] =",Integ[i])

Get_Integ = [0]*Close_count
for i in range(0,Close_count):
    Get_Integ[i] = integrate(Integ[i], (Xx, 0, Ll))

for i in range(0,Close_count):
    print("Get_Integ[",i+1,"] =",Get_Integ[i])

for i in range(0, len(Get_Integ)):
    print("I[",i,"] =",Get_Integ[i])
print("All I =",Get_Integ)
print("W_c =",W_coefs)
for solution in linsolve(Get_Integ, W_coefs):
    Result = solution

print(Result)

Result_num = []
for res in Result:
    Result_num.append(res.subs([(Ee,E), (Qq,q),(Ll,l),(Ii,I)]))

print("Coefs num =",Result_num)

Result_W_x = Get_Esp
Result_W_x = Result_W_x.subs([(Ll,l)])
for i in range(0,Close_count):
    Result_W_x = Result_W_x.subs([(W_coefs[i],Result_num[i])])

print("W(x) = ",Result_W_x)

x = np.linspace(0,l,Size_graph)
y_W = [0]*Size_graph

for i in range(0,Size_graph):
    y_W[i] = Result_W_x.subs([(Xx,x[i])])

plt.plot(x,y_W)
plt.show()

Result_M_x = -Ee*Ii*Dw[2]
Result_M_x = Result_M_x.subs([(Ll,l),(Ee,E),(Ii,I)])
for i in range(0,Close_count):
    Result_M_x = Result_M_x.subs([(W_coefs[i],Result_num[i])])

print("M(x) = ",Result_M_x)

y_M = [0]*Size_graph

for i in range(0,Size_graph):
    y_M[i] = Result_M_x.subs([(Xx,x[i])])

plt.plot(x,y_M)
plt.show()

Result_Q_x = Result_M_x.diff(Xx)

print("Q(x) = ", Result_Q_x)

y_Q = [0] * Size_graph

for i in range(0, Size_graph):
    y_Q[i] = Result_Q_x.subs([(Xx, x[i])])

plt.plot(x, y_Q)
plt.show()