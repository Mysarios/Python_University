import sympy
from sympy import *
from sympy import solveset, symbols, Interval, Min, Max
import numpy as np
import matplotlib.pyplot as plt
import random
import pylab
import scipy

Xx = Symbol('x')
def F_x_f(x):
    return 2 - x**2
def F_x_s(x):
    return -x

x_left = 0
x_right_first = 2**(1/2)
x_right_second = 2

y_down_f = 0
y_top_first = 2
y_top_second = -2

squid_f = x_right_first * y_top_first
squid_s = x_right_second * y_top_second * (-1)
m = 1000
Size = 50

good_f = 0
good_s = 0
bad_f = 0
bad_s = 0
X_now = - 5
X_array = []
lenght = 10
y_first = []
y_second = []
y_Line = []
for j in range(Size):
    y_Line.append(0)
    X_array.append(X_now)
    y_first.append(F_x_f(X_now))
    y_second.append(-F_x_s(X_now))
    X_now += lenght / Size
plt.plot(X_array, y_first)
plt.plot(X_array, y_second)
plt.show()
for i in range(m):
    x_f = random.random() * x_right_first
    x_s = random.random() * x_right_second

    y_f = random.random() * y_top_first
    y_s = random.random() * y_top_second

    Buf_y_f = F_x_f(x_f)
    Buf_s_f = F_x_s(x_s)
    Buf_s_f_2 = F_x_f(x_s)

    if Buf_y_f > y_f:
        good_f +=1
        plt.scatter(x_f, y_f, color='r', marker='o', s=20)
    else:
        bad_f +=1
        plt.scatter(x_f, y_f, color='b', marker='o', s=20)

    if Buf_s_f < y_s and Buf_s_f_2 > y_s:
        good_s +=1
        plt.scatter(x_s, y_s, color='r', marker='o', s=20)
    else:
        bad_s +=1
        plt.scatter(x_s, y_s, color='b', marker='o', s=20)

plt.plot(X_array, y_Line)
plt.show()

print(good_f/m)
print(good_s/m)

squid_first_fig = 2 * (squid_f * (good_f/m))
squid_second_fig =(squid_s * (good_s/m))
print(squid_second_fig)
print("First squire = ",squid_first_fig)
print("Second squire = ",squid_second_fig)


First_integral = 2 - Xx**2
First_integral = integrate(First_integral, (Xx, -x_right_first, x_right_first))
print("First int = ",First_integral)

Buf_inetgral = 2 - Xx**2
Buf_inetgral = integrate(Buf_inetgral, (Xx, -2, -x_right_first))

Second_integral = Xx
Second_integral = abs(integrate(Second_integral, (Xx, -2, 0)))
print("Second int = ",Second_integral - Buf_inetgral)

print("All squire = ",squid_first_fig + squid_second_fig)
print("All int = ",First_integral + Second_integral + Buf_inetgral)