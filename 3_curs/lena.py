import numpy as np
import matplotlib.pyplot as plt


y = [1]
x = np.linspace(0,19,20)
for i in range(0,len(x)-1):
    y.append(1 + 1/y[len(y)-1])

    print(round(y[i],7))
plt.plot(x,y)
plt.show()