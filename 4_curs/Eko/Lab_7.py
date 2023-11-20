import math
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
x_max = 120
y_max = 120
z_max = 540

X_mesh = []
Y_mesh = []
Z_mesh = []
for i in range(0,x_max+30,30):
    X_mesh.append(i)
    Y_mesh.append(i)
for i in range(0, z_max+30,30):
    Z_mesh.append(i)

start_t = 100
end_t = 2000
step_t = 10

high = 0
coordinates = 0
Stranght = 10
Disp_e1 = 1
Disp_e2 = 1
Disp_e3 = 1

# Мощность выброса
Stranght = 100
Mm = Stranght * 1000

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
        T = coef_turb_x/(S**2)
        return exp((-step_t)/(T))
    if i == 2:
        S = Var * i
        T = coef_turb_y/(S**2)
        return exp((-step_t)/(T))
    if i == 3:
        S = Var * i
        T = coef_turb_z/(S**2)
        return exp((-step_t)/(T))

def create_random(e1,d1,e2,d2,e3,d3,lenght):
    for i in range(0,int(lenght)):
        e1.append(random.random() * 2 * d1 - d1)
        e2.append(random.random() * 2 * d2 - d2)
        e3.append(random.random() * 2 * d3 - d3)


def create_club(e1,e2,e3,x,y,z,u,v,w,index,sigm_x1,sigm_x2,sigm_x3):
    u.append(middle_U + e1[index-1])
    v.append(middle_V + e2[index-1])
    w.append(middle_W + e3[index-1])
    x.append(0)
    y.append(0)
    z.append(30)

    sigm_x1.append(0)
    sigm_x2.append(0)
    sigm_x3.append(0)

def supply_step(x,y,z,u,v,w,index,t_now,sigm_x1,sigm_x2,sigm_x3):
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

    sigm_x1[index] = 2 * coef_turb_x * t_d
    sigm_x2[index] = 2 * coef_turb_y * t_d
    sigm_x3[index] = 2 * coef_turb_z * t_d

def get_concentration(x_arr,y_arr,z_arr,count_clubs,sigm_1,sigm_2,sigm_3,x,y):
    Concentrarion = 0
    for i in range(0,len(x)):
        #print("i =",i)
        for l in range(0,len(y)):
            #print("l =", l)
            for club in range(count_clubs):
                #print("club =", club)
                Concentrarion += (2 * Mm * step_t / ((2 * math.pi) ** 1.5 * np.sqrt(sigm_1[club] * sigm_2[club] * sigm_3[club]))
                       * math.exp(-(x[i] - x_arr[club]) ** 2 / (2 * sigm_1[club]))
                       * math.exp(-(y[l] - y_arr[club]) ** 2 / (2 * sigm_2[club]))
                       * math.exp(-z_arr[club] ** 2 / (2 * sigm_3[club])))

    return Concentrarion


def main_func():
    #start
    count_steps :int = (end_t-start_t) / step_t

    e_1 = []
    e_2 = []
    e_3 = []
    x_array = []
    y_array = []
    z_array = []
    u_array = []
    v_array = []
    w_array = []

    sigm_x1 = []
    sigm_x2 = []
    sigm_x3 = []

    x_for_graph = []
    y_for_graph = []


    result_concentration = []

    create_random(e_1,Disp_e1,e_2,Disp_e2,e_3,Disp_e3,count_steps)

    #включаем
    now_t = start_t
    count_club = 1
    create_club(e_1,e_2,e_3,x_array,y_array,z_array,u_array,v_array,w_array,count_club,sigm_x1,sigm_x2,sigm_x3)
    while now_t < end_t-step_t:
        print("Now_time",now_t)
        count_club += 1
        create_club(e_1, e_2, e_3, x_array, y_array, z_array, u_array, v_array, w_array,count_club,sigm_x1,sigm_x2,sigm_x3)
        #print("Create")
        for i in range(0,count_club):
            supply_step(x_array,y_array,z_array,u_array,v_array,w_array,i,now_t,sigm_x1,sigm_x2,sigm_x3)
            #print("Step")
        result_concentration.append(get_concentration(x_array,y_array,z_array,count_club,sigm_x1,sigm_x2,sigm_x3,X_mesh,Y_mesh))
        #print("Get_ conc")
        x_for_graph.append(now_t)
        y_for_graph.append(result_concentration[len(result_concentration) - 1])
        now_t+=step_t
    print(result_concentration)
    plt.plot(x_for_graph,y_for_graph)
    plt.show()

main_func()