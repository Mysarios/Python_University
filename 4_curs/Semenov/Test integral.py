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

quadratic_expression = 5*Xx + 5*Yy*Tt
#f = lambdify([Xx,Yy,Tt],[quadratic_expression])
f = lambda y, x,t: x + y*t
print(integrate.dblquad(f, 0, 1, 0, 1)[0])

quadratic_expression = 5*Xx + 5*Yy
quadratic_expression = factor_terms(quadratic_expression)
print(quadratic_expression)
sampling_points = np.linspace(low,high,n)
def inf_integrate(fun):

    f = lambdify([Xx,Yy],[fun])
    values = []
    for i in range(0,n):
        values.append(f(sampling_points[i],sampling_points[i]))
    print(values)
    return trapz(values,sampling_points)

#I1 = inf_integrate(quadratic_expression)
#print(I1)


