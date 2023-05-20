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
P_zero = []
P_zero.append([0.1,0.6,0.2,0,0.1])
Steps = 10

Matrix_next_step = [ [0.7,0.15,0.05,0,0.1],
                     [0.36,0.1,0.24,0.2,0.1],
                     [0,0.55,0.2,0.15,0.1],
                     [0.3,0.2,0.2,0.1,0.2],
                     [0.4,0.1,0.35,0.15,0] ]

for i in range(Steps):
    buf_step = np.matmul(P_zero[i-1],Matrix_next_step)
    P_zero.append(buf_step)
#print(P_zero)



Y_arrays = [[],[],[],[],[]]
#print(Y_arrays)
X_array = []
for i in range(Steps):
    Y_arrays[0].append(P_zero[i-1][0])
    Y_arrays[1].append(P_zero[i-1][1])
    Y_arrays[2].append(P_zero[i-1][2])
    Y_arrays[3].append(P_zero[i-1][3])
    Y_arrays[4].append(P_zero[i-1][4])
    X_array.append(i)

print(Y_arrays)
for i in range(5):
    plt.plot(X_array,Y_arrays[i-1])

plt.legend(["1","2","3","4","5"])
plt.show()
S_arr = []
Count_try = 5
for i in range(Count_try):
    S = 0
    for i in range(Steps):
        x=random.random()
        check = 1
        Sum = P_zero[i-1][S]

        while check:
            if x > Sum:

                if S !=4:
                    S += 1
                Sum += P_zero[i-1][S]
            else:
                check = 0
    S_arr.append(S)
    print(S)
y_array = []
count = 2
for i in range(Count_try):
    y_array.append(S_arr[0])
    for j in range(count):
        y_array[i-1] += S_arr[j]
    y_array[i-1] /=count
    print(y_array[i-1])

plt.plot(X_array,y_array)
plt.show()