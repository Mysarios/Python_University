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
Steps = 50

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

for i in range(5):
    plt.plot(X_array,Y_arrays[i-1])

plt.legend(["1","2","3","4","5"])
plt.show()

S_arr = []
X_arr = []
Count_try = 50

Counts = [0]*5
Count_S = 0
for i in range(Count_try):
    S = 1
    X_arr.append(i)
    Count_steps = 0
    while S!= 5:
        Count_steps += 1
        x=random.random()
        Sum = Matrix_next_step[S-1][0]
        check = 1
        ind = 1
        S_buf = 1
        while check:

            if x > Sum:
                S_buf += 1
                Sum += Matrix_next_step[S-1][ind]
                ind += 1
            else:
                check = 0
                S=S_buf
        Count_S+=1
        Counts[S-1] += 1
    S_arr.append(Count_steps)

print("P(1)= ",Counts[0]/Count_S)
print("P(2)= ",Counts[1]/Count_S)
print("P(3)= ",Counts[2]/Count_S)
print("P(4)= ",Counts[3]/Count_S)
print("P(5)= ",Counts[4]/Count_S)
print("Ter P(1)= ",Y_arrays[0][Steps-1])
print("Ter P(2)= ",Y_arrays[1][Steps-1])
print("Ter P(3)= ",Y_arrays[2][Steps-1])
print("Ter P(4)= ",Y_arrays[3][Steps-1])
print("Ter P(5)= ",Y_arrays[4][Steps-1])

y_array = []
count = 2
for i in range(Count_try):
    y_array.append(S_arr[0])
    for j in range(count):
        y_array[i-1] += S_arr[j-1]
    y_array[i-1] /= count
    count +=1


print(sum(S_arr)/Count_try)
print(y_array)
X_arr = X_arr[:-1]
y_array = y_array[:-1]
plt.plot(X_arr,y_array)
plt.show()

P1 = Symbol('p1')
P2 = Symbol('p2')
P3 = Symbol('p3')
P4 = Symbol('p4')
P5 = Symbol('p5')

P_array = [P1,P2,P3,P4,P5]
print(P_array)

P_next_array = [0]*6

count = 0
for j in range(5):
    for i in range(5):
        P_next_array[count] += Matrix_next_step[i][j]*P_array[i]

    count +=1

for j in range(5):
    P_next_array[j] -=P_array[j]
P_next_array[5] -= 1
for j in range(5):
    P_next_array[5] += P_array[j]


print(solve(P_next_array,P_array))



