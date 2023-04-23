import numpy
import numpy as np
import random
import matplotlib.pyplot as plt
from statistics import mean, variance
import math
import pylab


def f(x):
    return 2 / (math.pi * (1 - x**2)**(1/2))

Start = 0;
End = 1;

fig, ax = pylab.subplots()
ax.remove()

# print("Задание 1")
lbd = 2
m = 1000

pylab.subplot(2, 2, 1)
pylab.title("Плотность")
x = np.arange(Start, 5, 0.01)
plt.plot(x, lbd*math.e**(-lbd*x))

pylab.subplot(2, 2, 2)
pylab.title("Функция распределения")
x = np.arange(Start, 5, 0.01)
plt.plot(x, 1 - math.e**(-lbd*x))

data = []
pylab.subplot(2, 2, 3)
for i in range(m):
    data.append(-math.log(random.uniform(0, 1))/lbd)
    pylab.scatter(i + 1, data[i], color='r', marker='+', s=50)

pylab.subplot(2, 2, 4)
pylab.hist(data)

pylab.show()

print("Мат ожидание (теор):", 1/lbd)
print("Мат ожидание:", mean(data))
print("Дисперсия (теор):", 1/lbd**2)
print("Дисперсия:", variance(data))
print("Стандартное отклонение (теор):", 1/lbd)
print("Стандартное отклонение:", np.std(data))


print("Задание 2")
pylab.subplot(2, 2, 1)
pylab.title("Плотность")
x = np.arange(Start, End, 0.01)
plt.plot(x, 2 / (math.pi * (1 - x**2)**(1/2)))

pylab.subplot(2,2,2)
pylab.title("Функция распределения")
x = np.arange(0, 1, 0.01)
plt.plot(x, (2/math.pi) * numpy.arcsin(x))

data = []
pylab.subplot(2, 2, 3)
for i in range(m):
    data.append(random.uniform(0, 1)**(2/3))
    pylab.scatter(i + 1, data[i], color='r', marker='+', s=50)

pylab.subplot(2, 2, 4)
pylab.hist(data)

pylab.show()

print("Мат ожидание (теор):", 3*1**(5/2)/5)
print("Мат ожидание:", mean(data))
print("Дисперсия (теор):", 3*1**(7/2)/7 - (3*1**(5/2)/5)**2)
print("Дисперсия:", variance(data))
