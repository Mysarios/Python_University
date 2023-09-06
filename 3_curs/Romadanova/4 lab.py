import random

import sympy as sym
from sympy import *
import numpy as np
import math as m
import matplotlib as mpl
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

Data_X = [1,2,13]
Data_p = [0.1,0.1,0.8]

plt.title("Line graph")
plt.xlabel("X axis")
plt.ylabel("Y axis")
plt.plot(Data_X, Data_p, color ="red")
plt.show()

y_grapg_val = []
x_grapg_val = []
print("Ряд х = ",Data_X)
print("Ряд p этих х = ",Data_p)
print(" F(X<1) = 0")
print(" F(X<2) = 0.1")
print(" F(x<13) = 0.2")
print(" F(x<inf) = 1 ")
plt.title("Line graph")
plt.xlabel("X axis")
plt.ylabel("Y axis")

for i in range(-100,100):
    x_grapg_val.append(i)
    if i<=1:
        y_grapg_val.append(0)
    if i == 1:
        plt.plot(x_grapg_val, y_grapg_val, color ="red")
        y_grapg_val = []
        x_grapg_val = []
    if i <= 2 and i>1:
        y_grapg_val.append(0.1)
    if i == 2:
        plt.plot(x_grapg_val, y_grapg_val, color ="red")
        y_grapg_val = []
        x_grapg_val = []
    if i <= 13 and i>2:
        y_grapg_val.append(0.2)
    if i==13:
        plt.plot(x_grapg_val, y_grapg_val, color ="red")
        y_grapg_val = []
        x_grapg_val = []
    if i>13:
        y_grapg_val.append(1)
plt.plot(x_grapg_val, y_grapg_val, color ="red")
y_grapg_val = []
x_grapg_val = []

plt.plot(x_grapg_val, y_grapg_val, color ="red")
plt.show()

Count_rand = [0]*3
R_val =0
X_graph = []
Data_Analyz = []
for i in range(0,10000):
    X_graph.append(i)

    R_val = random.randint(1,10)
    if R_val == 1:
        Count_rand[0] +=1
        Data_Analyz.append(1)
    if R_val == 2:
        Count_rand[1] +=1
        Data_Analyz.append(2)
    if R_val > 2 :
        Count_rand[2]+=1
        Data_Analyz.append(13)

Count_rand[0]/=10000
Count_rand[1]/=10000
Count_rand[2]/=10000
print(Count_rand)
print(" Mean = ",np.array(Data_Analyz).mean())
print(" Var = ",np.array(Data_Analyz).var())
print(" Sdt = ",np.array(Data_Analyz).std())