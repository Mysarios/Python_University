import numpy as np
import random
import matplotlib.pyplot as plt
from statistics import mean, variance

import pylab

fig, ax = pylab.subplots()
ax.remove()

#Variant 15
Data_a = 5
Data_n = 12

Count_num=1000

DataBuf = []
for i in range(Count_num):
    DataBuf.append(random.uniform(Data_a, Data_a + Data_n ) // 1)

print(DataBuf)
pylab.subplot(2, 1, 1)

for i in range(Count_num):
    pylab.scatter(i + 1, DataBuf[i], color='b', marker='o', s=50)

pylab.subplot(2, 1, 2)
pylab.hist(DataBuf)

pylab.show()

Mean_Get = Data_a + (Data_n - 1) / 2
Variance_Get = (Data_n ** 2 - 1) / 12
print("Мат ожидание:")
print("По формуле:", Mean_Get)
print("По матлабу:", mean(DataBuf))
print("Дисперсия:")
print("По формуле:", Variance_Get)
print("По матлабу:", variance(DataBuf))
print("Стандартное отклонение:")
print("По формуле:", Variance_Get ** (1 / 2))
print("По матлабу:", np.std(DataBuf))
