from statistics import mean, variance
from sympy import *
from sympy import solveset, symbols, Interval, Min, Max
import numpy as np
import matplotlib.pyplot as plt
import random
import pylab


# g(x) = C*x**3  вариант 15 x in (0,2)
start = 0
end = 2
Xx = Symbol('x')
Cc = Symbol('С')

G_x = Cc*(Xx**3)
G_x_Integrated = integrate(G_x, (Xx, start, end))
C = 1 / (G_x_Integrated / Cc)

print("C = ",C)
print()

print("Find Max val of func")
lower_bound = start
upper_bound = end
f = C*Xx**3

zeros = solveset(f, Xx, domain=Interval(lower_bound, upper_bound))
assert zeros.is_FiniteSet # If there are infinite solutions the next line will hang.
res_min = Min(f.subs(Xx, lower_bound), f.subs(Xx, upper_bound), *[f.subs(Xx, i) for i in zeros])
res_max = Max(f.subs(Xx, lower_bound), f.subs(Xx, upper_bound), *[f.subs(Xx, i) for i in zeros])

print("res max = ",res_max)
print()

pylab.subplot(2, 2, 4)
x_array = np.arange(start, end, 0.01)
plt.plot(x_array, C*x_array**3)
#plt.show()


pylab.subplot(2, 2, 1)
print("Find eps and nj")
m=1000
data = []
data_for_val = []
x_array = np.arange(start, end, 0.01)
Max =0
for i in range(len(x_array)):
    pylab.scatter(x_array[i],C*x_array[i]**3, color='g', marker='o', s=20)
for i in range(m):
    x1 = random.uniform(start, end)
    x2 = random.uniform(start, end)

    y1 = C*x1**3
    y2 = C*x2**3

    eps =x1
    nj = x2
    gx = C*eps**3
    if nj> Max:
        Max = nj
    if eps <= 2:
        if nj >= gx:
            pylab.scatter(eps,nj, color='r', marker='+', s=50)
        else:
            pylab.scatter(eps, nj, color='b', marker='x', s=50)
            data.append(eps)
            data_for_val.append(eps)




#pylab.show()
pylab.subplot(2, 2, 2)
print("Max = ",Max)
data = np.array(data)

data_graph = []
data_y = []
for i in range(len(data)):
    pylab.scatter(i, data[i], color='r', marker='o', s=20)
    data_y.append(data_for_val[i])

data_graph = data.astype(np.float64)
data_y = np.array(data_y)
data_y = data_y.astype(np.float64)

#pylab.show()
pylab.subplot(2, 2, 3)
pylab.hist(data_graph)
pylab.show()

G_x = C*(Xx**3)*Xx
G_x_Integrated = integrate(G_x, (Xx, start, end))
Mean_Ter = G_x_Integrated

G_x = C*(Xx**3)*Xx*Xx
G_x_Integrated = integrate(G_x, (Xx, start, end))
Var_Ter = G_x_Integrated - Mean_Ter**2


print("Мат ожидание (теор):", Mean_Ter.evalf())
print("Мат ожидание:", mean(data_y))
print("Дисперсия (теор):", Var_Ter.evalf())
print("Дисперсия:", variance(data_y))
