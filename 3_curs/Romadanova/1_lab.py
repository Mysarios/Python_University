import numpy as np
import matplotlib.pyplot as plt
import random
from matplotlib import style
import statistics as st
import scipy


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
    a=5
    b=12
    c = []
    j = 0
    for i in range(len(x)):
        c.append(a + (b-a)*random.random())
        j += 1
    return c

plt.subplot(2,2,1)
x=np.arange(-22,22,0.01)
plt.plot(x,First_plot(x))
#plt.show()

plt.subplot(2,2,2)
plt.plot(x,Second_plot(x))
#plt.show()

plt.subplot(2,2,3)
style.use('ggplot')
X_third = []
for i in range(1 , 1000):
    X_third.append(i,)

g=Third_plot(X_third)
#print(g)
a=5
b=12
numbers= np.arange(a,b,(b-a)/10)
plt.scatter(X_third,g)
#plt.show()

plt.subplot(2,2,4)
plt.hist(g,numbers)
plt.show()

a,b=0,1
print("Math E = " + str(st.mean(g)) + " "  + str((a+b)/2))
print("Math D = " + str(np.var(g))  + " "  + str(((b-a)**2)/12))
print("Math D^1/2 = " + str(st.stdev(g))+ " "+ str((((b-a)**2)/12)**(1/2.0)))


print("a = ", end=" ")
a = int(input())
print("b =", end=" ")
b = int(input())
print("m = ", end=" ")
m = int(input())

data = []
for i in range(m):
    data.append(random.uniform(a, b))

print()
print("Question 1:")
#print(data)
print()
print("Question 2:")
x = random.uniform(a, b)
print(x, scipy.stats.uniform.cdf(x, a, b - a))
print()
print("Question 3:")
x1 = random.uniform(a, b)
x2 = random.uniform(a, b)
print(x1, x2, abs(scipy.stats.uniform.cdf(x1, a, b - a) - scipy.stats.uniform.cdf(x2, a, b - a)))
print()
print("Question 4:")
print(scipy.stats.uniform.pdf(x, a, b - a), scipy.stats.uniform.pdf(a - 1, a, b - a))
print()
print("Question 5:")
print(scipy.stats.uniform.ppf(0.5, a, b - a))
# print(scipy.stats.norm.ppf(0.5,a,b-a))
print()
print("Question 6:")
print(scipy.stats.uniform.stats(a, b - a))
print(scipy.stats.uniform.mean(a, b - a))
print(scipy.stats.uniform.var(a, b - a))


