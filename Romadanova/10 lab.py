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


m = 1000

n = 60
X = [0]*m
R = [0]*n
X_array = []

mu = 8.6
sigma = 2.4
for i in range(m):
    X_array.append(i)
    #s = np.random.normal(mu, sigma, n)
    for j in range(n):
        R[j] = random.random()
    X[i] = (sum(R) - (n/2))/((n/12)**(1/2))*(sigma**(1/2)) + mu

pylab.subplot(1, 2, 1)
plt.scatter(X_array,X)
pylab.subplot(1, 2, 2)
pylab.hist(X)
plt.show()
print("Question 1")
print("Мат ожидание (теор):", mu)
print("Мат ожидание:", mean(X))
print("Дисперсия (теор):", sigma)
print("Дисперсия:", variance(X))

mu = 8.6
sigma = 2.4

r =[0]*m
fi = [0]*m
z0 = [0]*m
z1 = [0]*m

for j in range(m):
    r = random.random()
    fi = random.random()
    if r == 0:
        r+= 0.0000000001
    if fi == 0:
        fi += 0.0000000001

    z0[j] = cos(2*math.pi*fi) * ((-2*math.log10(r))**(1/2))
    z1[j] = sin(2*math.pi*fi) * ((-2*math.log10(r))**(1/2))
    z0[j] = z0[j]*(sigma**(2/2)) + mu
    z1[j] = z1[j] * (sigma ** (2 / 2)) + mu


Data_z0 = z0.copy()
Data_z1 = z1.copy()
pylab.subplot(2, 2, 1)
plt.scatter(X_array,z0)
pylab.subplot(2, 2, 2)
z0 = np.array(z0)
z0 = z0.astype(np.float64)
pylab.hist(z0)
pylab.subplot(2, 2, 3)
z1 = np.array(z1)
z1 = z1.astype(np.float64)
plt.scatter(X_array,z1)
pylab.subplot(2, 2, 4)
pylab.hist(z1)
plt.show()

print("Question 2")
print("Мат ожидание (теор):", mu)
print("Мат ожидание:", mean(z0))
print("Дисперсия (теор):", sigma)
print("Дисперсия:", variance(z0))

print("Question 3")
print("Мат ожидание (теор):", mu)
print("Мат ожидание:", mean(z1))
print("Дисперсия (теор):", sigma)
print("Дисперсия:", variance(z1))
