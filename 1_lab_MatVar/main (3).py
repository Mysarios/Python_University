import numpy as np
import matplotlib.pyplot as plt
from random import randint
from statistics import mean

import pylab

fig, ax = pylab.subplots()
ax.remove()

n = 1000
m = [1291]
M = 21870
c = 0
g = 12
gamma = [m[0] / M]
period = []

pylab.subplot(1, 2, 1)
for i in range(n):
    m.append((g * m[i] + c) % M)
    gamma.append(m[i + 1] / M)
    pylab.scatter(i + 1, gamma[i + 1], color='r', marker='+', s=50)
    if gamma[i] == gamma[1]:
        period.append(i - 1)

gamma.pop(0)

pylab.subplot(1, 2, 2)
pylab.hist(gamma, bins=10)
pylab.suptitle("c = 0")
pylab.show()
print("c =", c)
print("Per = ", period[1])
print("n =", n)
print("Mat_Exp M(X) = %.5f" % mean(gamma))
print("D D(X) = %.5f" % np.var(gamma))
print()

n = 1000
m = [1291]
M = 21870
g = 12
c = 4621
gamma = [m[0] / M]

pylab.subplot(1, 2, 1)
for i in range(n):
    m.append((g * m[i] + c) % M)
    gamma.append(m[i + 1] / M)
    pylab.scatter(i + 1, gamma[i + 1], color='r', marker='+', s=50)
    if gamma[i] == gamma[1]:
        period.append(i - 1)

gamma.pop(0)

pylab.subplot(1, 2, 2)
pylab.hist(gamma, bins=10)
pylab.suptitle("c = %i" % c)
pylab.show()

print("c =", c)
print("Per = ", period[1])
print("n =", n)
print("Mat_Exp M(X) = %.5f" % mean(gamma))
print("D D(X) = %.5f" % np.var(gamma))
print()