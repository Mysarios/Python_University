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


# g(x) = C*x**3  вариант 15 x in (0,2)
start = -3
end = 3
Xx = Symbol('x')
Tt = Symbol('t')



nu = 0
sigma =2
m = 1000

f_x = (1 /(2*math.pi*sigma**2)**1/2) * exp( - ((Xx - nu)**2)/ 2*sigma**2)
x_array = np.arange(start, end, 0.01)
y_array = []
for i in range(len(x_array)):
    y_array.append(f_x.subs(Xx,x_array[i]))

pylab.subplot(2, 2, 1)
plt.plot(x_array,y_array)

F_x = (1 /(2*math.pi*sigma**2)**1/2) * exp( - ((Tt - nu)**2)/ 2*sigma**2)
F_x = integrate(F_x, (Tt, 0, Xx))

y_array_2 = []
for i in range(len(x_array)):
    y_array_2.append(F_x.subs(Xx,x_array[i]))


pylab.subplot(2, 2, 2)
plt.plot(x_array,y_array_2)

#r = norm.rvs(size=m)
r = np.random.normal(nu, sigma, m)
pylab.subplot(2, 2, 3)
x_data = []
y_data = []
for i in range(m):
    x = random.uniform(start, end)
    y = F_x.subs(Xx,x)
    x_data.append(i)
    y_data.append(y)
    #pylab.scatter(i,y, color='b', marker='o', s=6)

pylab.scatter(x_data,r, color='b', marker='o', s=6)


pylab.subplot(2, 2, 4)
pylab.hist(r)
plt.show()




print("Мат ожидание (теор):", nu)
print("Мат ожидание:", mean(r))
print("Дисперсия (теор):", sigma**2)
print("Дисперсия:", variance(r))