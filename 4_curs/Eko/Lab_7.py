from numba import njit
from sympy import *
import sympy as sym
import numpy as np
import math as m
import matplotlib.pyplot as plt
import random


#Входные

#Symbols
Xx = Symbol('x')
Yy = Symbol('y')
Zz = Symbol('z')

x_point = 0
y_point = 0
z_point = 30

i=0.1

x_max = 9000
y_max = 9000
z_max = 500

start_t = 100
end_t = 3700
step_t = 10

high = 0
coordinates = 0
Stranght = 10
Disp_e1 = 1
Disp_e2 = 1
Disp_e3 = 1

middle_U = 3
middle_W = 0.3
middle_V = 3

coef_turb_x = 100
coef_turb_y = 1000
coef_turb_z = 10
# functions
def get_Pul(i,Var):
    if i == 1:
        S = Var * i
        T = coef_turb_x/(S^2)
        exp((-step_t)/(T))
    if i == 2:
        S = Var * i
        T = coef_turb_y/(S^2)
        exp((-step_t)/(T))
    if i == 3:
        S = Var * i
        T = coef_turb_z/(S^2)
        exp((-step_t)/(T))

def create_random(e1,d1,e2,d2,e3,d3,lenght):
    for i in range(0,int(lenght)):
        e1.append(random.random() * 2 * d1 - d1)
        e2.append(random.random() * 2 * d2 - d2)
        e3.append(random.random() * 2 * d3 - d3)
def create_club(e1,e2,e3,x,y,z,u,v,w):
    u.append(middle_U + e1)
    v.append(middle_V + e2)
    w.append(middle_W + e3)
    x.append(0)
    y.append(0)
    z.append(30)
def supply_step(x,y,z,u,v,w,index,t_now):
    P_1 = get_Pul(1, u[index])
    P_2 = get_Pul(2, v[index])
    P_3 = get_Pul(3, w[index])

    u[index] = u[index] * P_1 + random.random() * 2 * Disp_e1 - Disp_e1
    v[index] = v[index] * P_2 + random.random() * 2 * Disp_e2 - Disp_e2
    w[index] = w[index] * P_3 + random.random() * 2 * Disp_e3 - Disp_e3

    x[index] += u[index]*step_t
    y[index] += v[index]*step_t
    z[index] += w[index]*step_t
    t_d = t_now - step_t*index

    Sigma_x1 = 2 * coef_turb_x * t_d
    Sigma_x2 = 2 * coef_turb_y * t_d
    Sigma_x3 = 2 * coef_turb_z * t_d


#start
count_steps :int = (end_t-start_t) / step_t

E_1 = []
E_2 = []
E_3 = []
X_array = []
Y_array = []
Z_array = []
U_array = []
V_array = []
W_array = []

create_random(E_1,Disp_e1,E_2,Disp_e2,E_3,Disp_e3,count_steps)

#включаем
now_t = start_t
create_club(E_1,E_2,E_3,X_array,Y_array,Z_array,U_array,V_array,W_array)
count_club = 1
while now_t < end_t:
    create_club(E_1, E_2, E_3, X_array, Y_array, Z_array, U_array, V_array, W_array)
    count_club += 1
    for i in range(0,count_club):
        supply_step(X_array,Y_array,Z_array,U_array,V_array,W_array,i,now_t)

