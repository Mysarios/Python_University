import matplotlib.pyplot as plt
import numpy as np
from statistics import mean
import pylab


def f(x,K):
    Size=1
    print(x)
    buf = x*(10**(K/4))
    print(buf)
    buf-=(buf // 1)
    print(buf)
    buf *= 10 ** (K / 2)
    print(buf)
    buf = buf // 1
    print(buf)



    buf = 0 + (buf**2)*(10**(-1*K))
    print(buf)
    result = buf
    return result
def a_f(x,K):
    Size = 1
    buffer = x**2
    if (x*10000)**2 <= 10**15:
        buffer = buffer/10
    if (x*10000)**2 >= 10**16:
        buffer = buffer*10
    buffer = buffer*(10**(K/4))
    buffer-=(buffer // 1)
    buffer *= 10 ** (K / 2)
    buffer = buffer // 1
    buffer = 0 + (buffer)/((10**(K/2))*Size)
    result = buffer
    return result

n =1000
k=16
x_array = [0]
y_array_1 = [0.7234567899876543]
y_array = [0.56788765]
buf_num =0

for i in range(1,n):
    x_array.append(i)
    y_array.append(a_f(y_array[i - 1], k))
    if a_f(y_array[i-1],k) == 0:
        break


fig, ax = pylab.subplots()
ax.remove()
pylab.subplot(1, 2, 1)
plt.scatter(x_array,y_array)

pylab.subplot(1, 2, 2)
Y_range = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]
plt.hist(y_array,Y_range)
plt.show()


fig, ax = pylab.subplots()
ax.remove()

x_array = [0]
for i in range(1,n):
    x_array.append(i)

n = 1000
m = [1291]
M = 21870
c = 0
g = 16
gamma = [m[0] / M]
period = []

pylab.subplot(1, 2, 1)
for i in range(0,n-1):
    m.append((g * m[i] + c) % M)
    gamma.append(m[i + 1] / M)
    if gamma[i] == gamma[1]:
        period.append(i - 1)


pylab.scatter(x_array, gamma, color='g', marker='+')
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

m_2 = [1291]
g = 11
c = 4621
gamma = [m_2[0] / M]

pylab.subplot(1, 2, 1)
for i in range(0,n-1):
    m_2.append((g * m_2[i] + c) % M)
    gamma.append(m_2[i + 1] / M)
    pylab.scatter(i + 1, gamma[i + 1], color='g', marker='+')
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