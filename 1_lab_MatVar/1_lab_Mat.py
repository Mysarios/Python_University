import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import style


def First_plot(x):
    a=5
    b=12
    c = []
    j=0
    for i in range(len(x)):
        if (x[j]>a and x[j]<b):
            c.append(1/(b-a))
        else:
            c.append(0)
        j+=1
    return c
def Second_plot(x):
    a = 5
    b = 12
    c = []
    j = 0
    for i in range(len(x)):
        if (x[j] < a):
            c.append(0)
        elif(x[j]>=a and x[j]<b):
            c.append((x[j]-a)/(b-a))
        elif(x[j]>=b):
            c.append(1)

        j += 1
    return c
def Third_plot(x):
    a=0
    b=10
    c = []
    j = 0
    for i in range(len(x)):
        c.append(random.random())
        j += 1
    return c

x=np.arange(-22,22,0.01)
plt.plot(x,First_plot(x))
plt.show()


plt.plot(x,Second_plot(x))
plt.show()

style.use('ggplot')
X_third = []
for i in    range(1,1000):
    X_third.append(random.random())

g=Third_plot(X_third)
print(g)
numbers =[2 , 2.5 , 3 ,3.5, 4 , 4.5 , 5]
plt.scatter(X_third,g)
plt.show()

plt.hist(g,numbers)
plt.xlabel('percentage')
plt.ylabel('Number of people')
plt.title('Histogram')
plt.show()



