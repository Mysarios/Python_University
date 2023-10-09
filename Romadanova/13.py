import math
from statistics import mean, variance
from scipy. stats import norm
import sympy
from sympy import *
from sympy import solveset, symbols, Interval, Min, Max
import numpy as np
import matplotlib.pyplot as plt
import random
import pylab
import scipy


Steps = 10

Matrix_next_step = [ [0,-1,-1/3,0],
                     [-1/2,0,0,-1],
                     [-1/3,0,0,-1],
                     [0,-1/2,-1/2,0]]

S = 1
x_array = []
y_array = []
LNR = [0,0,0,0]
y_array.append(S)
x_array.append(0)
for i in range(Steps):
    x_1 = random.random()
    x_2 = random.random()
    x_3 = random.random()
    x_4 = random.random()
    LNR = [0, 0, 0, 0]
    LNR[0] = math.log(x_1)
    LNR[1] = math.log(x_2)
    LNR[2] = math.log(x_3)
    LNR[3] = math.log(x_4)
    min = 100
    S_buf = S
    for j in range(4):
        if Matrix_next_step[S-1][j] * LNR[j] < min and Matrix_next_step[S-1][j] * LNR[j] != 0:
            min = Matrix_next_step[S - 1][j] * LNR[j]
            S_buf = j+1
            print(f"Step = {i}", f"S = {S_buf}", f"Min val = {min}")

    S = S_buf
    y_array.append(S)
    x_array.append(i)

print(y_array)
print(x_array)
plt.plot(x_array,y_array)
plt.show()
