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
x_max = 9000
y_max = 9000
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
end_t = 1000
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
middle_W = 0.003
middle_V = 0.0003

coef_turb_x = 100
coef_turb_y = 100
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
        e1 = np.append(e1,random.random() * 2 * d1 - d1)
        e2 = np.append(e2,random.random() * 2 * d2 - d2)
        e3 = np.append(e3,random.random() * 2 * d3 - d3)
    return e1,e2,e3

def create_club(e1,e2,e3,x,y,z,u,v,w,index,sigm_x1,sigm_x2,sigm_x3):
    u = np.append(u,middle_U + e1)
    v = np.append(v,middle_V + e2)
    w = np.append(w,middle_W + e3)
    x = np.append(x,0)
    y = np.append(y,0)
    z = np.append(z,30)

    sigm_x1 = np.append(sigm_x1,0)
    sigm_x2 = np.append(sigm_x2,0)
    sigm_x3 = np.append(sigm_x3,0)

    return u,v,w,x,y,z,sigm_x1,sigm_x2,sigm_x3

def supply_step(x,y,z,u,v,w,index,t_now,sigm_x1,sigm_x2,sigm_x3):
    P_1 = get_Pul(1, u[index])
    P_2 = get_Pul(2, v[index])
    P_3 = get_Pul(3, w[index])


    #u[index] = u[index] * P_1 + random.random() * 2 * Disp_e1 - Disp_e1
    #v[index] = v[index] * P_2 + random.random() * 2 * Disp_e2 - Disp_e2
    #w[index] = w[index] * P_3 + random.random() * 2 * Disp_e3 - Disp_e3

    u[index] = u[index] * P_1 + random.random()
    v[index] = v[index] * P_2 + random.random()
    w[index] = w[index] * P_3 + random.random()

    x[index] += u[index]*step_t
    #print("New x =",x[index],"  u = ",u[index])
    y[index] += v[index]*step_t*0.0001
    z[index] += w[index]*step_t
    t_d = t_now - step_t*index

    sigm_x1[index] = (2 * coef_turb_x * t_d)**(1/2)
    sigm_x2[index] = (2 * coef_turb_y * t_d)**(1/2)
    sigm_x3[index] = (2 * coef_turb_z * t_d)**(1/2)

    return u,v,w,x,y,z,sigm_x1,sigm_x2,sigm_x3
@njit
def get_concentration(x_arr,y_arr,z_arr,count_clubs,sigm_1,sigm_2,sigm_3,x,y):
    Concentrarion = 0
    for i in range(0,len(x)):
        #print("i =",i)
        for l in range(0,len(y)):
            #print("l =", l)
            for club in range(count_clubs):
                #print("club =", club)
                Concentrarion += (2 * Mm * step_t / ((2 * math.pi) ** 1.5 * np.sqrt(sigm_1[club] * sigm_2[club] * sigm_3[club])) * math.exp(-(x[i] - x_arr[club]) ** 2 / (2 * sigm_1[club]))* math.exp(-(y[l] - y_arr[club]) ** 2 / (2 * sigm_2[club])))

    return Concentrarion
@njit
def get_concentration_x(Buffer,x_arr,y_arr,z_arr,count_clubs,sigm_1,sigm_2,sigm_3,x,y):

    #print("i =",i)
    for ind in range(0,300):
        Buffer[ind] = Buffer[ind]*(count_clubs-1)
    for i in range(0,len(x)):
        for j in range(0,len(y)):
            for club in range(count_clubs):
                #print("Club =",club,"  x=",x_arr[club])
                #print("club =", club)
                Buffer[i] += (Mm * step_t / (2 * math.pi * sigm_1[club] * sigm_2[club]))\
                        * math.exp(-((x[i] - x_arr[club]) ** 2) / (2 * sigm_1[club]**2) - ((y[j] - (y_arr[club])) ** 2) / (2 * sigm_2[club]**2)) \
                             * (2/(( (2*math.pi)**(1/2))*sigm_3[club])) * math.exp(((- 30) ** 2) / (2 * sigm_3[club]**2))
    for ind in range(0, 300):
        Buffer[ind] /= count_clubs
                       #* math.exp(( 0** 2)/ (2 * sigm_3[club])))
    #print("x cood by 1 club",x_arr[0],"  y coord by 1 club =",y_arr[0])
    #print("conc =",0,"  = ", Buffer[0])
    #print("conc =", 10, " = ", Buffer[10])
    #print("conc =", 20, " = ", Buffer[20])

    #print("Conc =",Buffer)
    return Buffer


def main_func():
    #start
    count_steps :int = (end_t-start_t) / step_t

    Concentrarion_x = np.array([])
    x_array = np.array([])
    y_array = np.array([])
    z_array = np.array([])
    u_array = np.array([])
    v_array = np.array([])
    w_array = np.array([])

    sigm_x1 = np.array([])
    sigm_x2 = np.array([])
    sigm_x3 = np.array([])

    x_for_graph = np.array([])
    y_for_graph = np.array([])
    x_for_graph_2 = np.array([])
    y_for_graph_2 = np.array([])


    result_concentration = np.array([])
    result_concentration_x = np.array([])
    for i in range(0,300):
    #array_results_x = np.array([])
        result_concentration_x = np.append(result_concentration_x,0)

    e_1 = random.random() * 2 * Disp_e1 - Disp_e1
    e_2 = random.random() * 2 * Disp_e2 - Disp_e2
    e_3 = random.random() * 2 * Disp_e3 - Disp_e3
    #включаем
    now_t = start_t
    count_club = 1
    u_array,v_array,w_array,x_array,y_array,z_array,sigm_x1,sigm_x2,sigm_x3 = create_club(e_1,e_2,e_3,x_array,y_array,z_array,u_array,v_array,w_array,count_club,sigm_x1,sigm_x2,sigm_x3)
    result_concentration_x = np.append(result_concentration_x, 0)
    while now_t < end_t-step_t:
        print("Now_time",now_t)
        count_club += 1
        e_1 = random.random() * 2 * Disp_e1 - Disp_e1
        e_2 = random.random() * 2 * Disp_e2 - Disp_e2
        e_3 = random.random() * 2 * Disp_e3 - Disp_e3
        u_array,v_array,w_array,x_array,y_array,z_array,sigm_x1,sigm_x2,sigm_x3 = create_club(e_1, e_2, e_3, x_array, y_array, z_array, u_array, v_array, w_array,count_club,sigm_x1,sigm_x2,sigm_x3)
        for i in range(0,count_club):
            u_array,v_array,w_array,x_array,y_array,z_array,sigm_x1,sigm_x2,sigm_x3 = supply_step(x_array,y_array,z_array,u_array,v_array,w_array,i,now_t,sigm_x1,sigm_x2,sigm_x3)
        result_concentration = np.append(result_concentration,get_concentration(x_array,y_array,z_array,count_club,sigm_x1,sigm_x2,sigm_x3,X_mesh,Y_mesh))
        result_concentration_x = get_concentration_x(result_concentration_x,x_array, y_array, z_array, count_club, sigm_x1, sigm_x2,sigm_x3, X_mesh, Y_mesh)
        x_for_graph = np.append(x_for_graph,now_t)
        y_for_graph = np.append(y_for_graph,result_concentration[len(result_concentration) - 1])
        now_t+=step_t

    result_x = 0
    for i in range(0,x_max+30,30):
        x_for_graph_2 = np.append(x_for_graph_2, i)

    plt.plot(x_for_graph,y_for_graph)
    plt.show()
    plt.plot(x_for_graph_2, result_concentration_x)
    plt.show()

    xm = 5
    print(xm)

    Cm = (Mm/(30*30))
main_func()