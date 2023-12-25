from threading import Thread
from symengine import *
import sympy
from scipy.integrate import trapz
import numpy as np
from sympy import factor_terms
from scipy import integrate

Xx, Yy = symbols('x, y')
Tt = symbols('t')
low = 0
high = 1
n = 5000
fun = []
quadratic_expression_1 = 5*Xx + 5*Yy
quadratic_expression_2 = 22*Xx**10 + 5*Yy
fun.append(quadratic_expression_1)
fun.append(quadratic_expression_2)

def intagrete(fun,i):
    fun[i] = integrate(fun[i], (Xx, 0, 1))

print(fun[0])
print(fun[1])

t1 =Thread(target=intagrete, args=(fun,0,))
t2 =Thread(target=intagrete, args=(fun,1,))

t1.start()
t2.start()

t1.join()
t2.join()

print(fun[0])
print(fun[1])
