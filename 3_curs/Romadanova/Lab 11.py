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


def f(x):
    return (x**3)*((0.7 + x)**(1/2))
def f_2(x,y):
    return 6*x*(y**2)-4*x
n = 10000

Xx = Symbol('x')
Yy = Symbol('y')

a = 1
b = 2.2
Disp = b-a

eps = [0]*n
f_eps = [0]*n
for i in range(n):
    eps[i] = random.random()*Disp + a
    f_eps[i] = f(eps[i])
Integral_val_MK = ((b-a)/n)*sum(f_eps)
print("Our num =",Integral_val_MK)
print("teor num = ", scipy.integrate.quad(f, a, b)[0])



a_x = 1
b_x = 2
a_y = 2
b_y = 3


n = 10000

Xx = Symbol('x')
Yy = Symbol('y')

Disp_ab = b_x-a_x
Disp_cd = b_y-a_y

x = [0]*n
y = [0]*n
f_eps = [0]*n

F_x = 6*Xx*(Yy**2)-4*Xx
F_x = integrate(F_x, (Xx, a_x, b_x))
F_x = integrate(F_x, (Yy, a_y, b_y))

for i in range(n):
    x[i] = random.random() * Disp_ab + a_x
    y[i] = random.random() * Disp_cd + a_y
    f_eps[i] = f_2(x[i],y[i])
Integral_val_MK = ((b_x-a_x)*(b_y-a_y)/n)*sum(f_eps)
print("Our num =",Integral_val_MK)
print("teor num = ", F_x)