import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
from numpy import linalg as LA
from sympy import *

K = 8
L = 3
size = K+1
u = 5
l = 2.5

Lyambda = l
Lyambda_2 = u

Array_intens = [0]*(size)
for i in range (0,len(Array_intens)):
    Array_intens[i] = [0]*(size)


#Uu = Symbol('u')
intens = 1
for low in range(1,len(Array_intens)):
    Array_intens[low][low-1] = intens * u
    if intens < L:
        intens +=1

#Ll = Symbol('l')
for high in range(0,len(Array_intens)-1):
    Array_intens[high][high+1] = l*(high+1)


for middle in range(0,len(Array_intens)-1):
    Array_intens[middle][middle] -= (Array_intens[middle][middle-1] + Array_intens[middle][middle+1])

Array_intens[len(Array_intens)-1][len(Array_intens)-1] = -L*u

print("base =")
for i in range(0,size):
    print(Array_intens[i])

Matrix = np.array(Array_intens)
Matrix = Matrix.transpose()
print("transp base = =")
for i in range(0,size):
    print(Matrix[i])
eigenvalues, eigenvectors = LA.eig(Matrix)
print("E values = ")

for i in range(0,size):
    if abs(eigenvalues[i]) - 0.005 < 0:
        eigenvalues[i] = 0
    print("E" + str(i) + " ",eigenvalues[i])

print("E vectors = ")

for i in range(0, size):
    buf = "EV"+ str(i)
    for j in range(0, size):
        buf += ": {:.4f}".format(eigenvectors[i][j]) + " "
    print(buf)


C_symbols =[0]*(size)
P_symbols = [0]*(size)
arrayFor_K_CH = [0]*(size)
Static_eval = [0]*(size)
for i in range (0,size):
    #arrayFor_K_CH[i] = [0]*9
    #Static_eval[i] = [0]*9
    C_symbols[i] = Symbol('C' + str(i))
    P_symbols[i] = Symbol('P' + str(i))

arrayFor_K_CH[0] = 0
Static_eval[0] = -1

Tt = Symbol('t')
for i in range(0, size):
    for j in range(0,size):
        #arrayFor_K_CH[i] += C_symbols[j] * eigenvectors[i][j] * math.e**(Tt * eigenvalues[j])
        Static_eval[i] += C_symbols[j] * eigenvectors[i][j]

print("System with C  = ")
for i in range(0, size):
    buf = 0
    for j in range(0, size):
        buf += C_symbols[j] * round(eigenvectors[i][j],2) * Symbol('e')**(Tt * round(eigenvalues[j],1))
    print(buf)


Result = 0
for solution in linsolve(Static_eval, C_symbols):
    Result = solution

print("C coefs  =",Result)


Transp_Array_intens = Array_intens.copy()
Transp_Array_intens = np.array(Transp_Array_intens)
Transp_Array_intens = Transp_Array_intens.transpose()

for i in range(0,size):
    print(Transp_Array_intens[i])

print("System   = ")
for i in range(0, size):
    buf = 0
    for j in range(0, size):
        arrayFor_K_CH[i] += Result[j] * eigenvectors[i][j] * math.e**(Tt * eigenvalues[j])
        buf += round(Result[j] * eigenvectors[i][j],4) * Symbol("e")**(Tt * round(eigenvalues[j],1))
    print(buf)

print("Results   = ")
Sym = 0
limit_expr = []
for i in range (0,len(arrayFor_K_CH)):
    limit_expr.append(limit(arrayFor_K_CH[i], Tt, oo))
    Sym += limit_expr[i]
    print(" Limit ",i," = ",round(limit_expr[i],9))
print("Check sym =", round(Sym,2))


print("Effectivnes params by t = oo:")
systemStrain = Lyambda/Lyambda_2
print(" Strain = ",systemStrain)
systemStrainPerChannel = systemStrain/2
print(" Strain per channel = ",systemStrainPerChannel)

chanceOfFail = 0
print(" Chance fail = ", chanceOfFail)
chanceBuzy = 1 - round(limit_expr[0],5)
print(" Chance Service = ", chanceBuzy)
relativeThroughput = 1
print(" Relative throughput = ", relativeThroughput)
absoluteThroughput = Lyambda_2 * chanceBuzy
print(" Absolute throughput = ", absoluteThroughput)

averageMaintenance = 0                              #Среднее заявок под обслуживанием  Nоб
for k in range(0,L):
    averageMaintenance += k * round(limit_expr[k],5)
for k in range(L,size):
    averageMaintenance += L * round(limit_expr[k],5)
print(" ")
print(" Average maintenance = ", averageMaintenance)

averageSystem = 0                                   #Среднее пасивов  Nсист
for k in range(0,size):
    averageSystem += k * round(limit_expr[k],5)
print(" Average system = ", averageSystem)

averageQueue = averageSystem - averageMaintenance   #Среднее заявок  в очереди Nоч
print(" Average queue = ", averageQueue)

averageWait = 0                                     #Среднее простоя Nпр
for k in range(0,L):
    averageWait += (L-k) * round(limit_expr[k],5)
print(" Average wait = ", averageWait)

Delta = (K - averageSystem)* Lyambda  # Дельта

averageWaitQueue_Time = (1/Delta) * averageQueue    #Среднее время в очереди Tоч
print(" Average wait Queue Time = ", averageWaitQueue_Time)

averageServiceTime = (1 / Delta) * averageMaintenance                    #Среднее время обслуживания Tоб
print(" Average Service Time = ", averageServiceTime)

averageSystemTime = averageServiceTime + averageWaitQueue_Time  #Среднее время в системе Tсист
print(" Average System Time = ", averageSystemTime)

print("Formuls Littl:")

T_Maintenance = (1 / Delta) * averageQueue #Среднее время в очереди Tоч
print(" T_wait = ", T_Maintenance)

T_Servise = (1 / Delta) * averageMaintenance   #Среднее время обслуживания Tоб
print(" T_Service Time = ", T_Servise) #&&

T_System = averageSystem/Delta               #Среднее время в системе Tсист
print(" T_System = ", T_System)




print("Effectivnes params by dynamic t :")
time = 20*5
timeForPlot = []
variantsPerTime = [0] * time
for t in range (0,time):
    timeForPlot.append(t/5)
    variantsPerTime[t] = [0] * len(arrayFor_K_CH)
    for p in range(0,len(arrayFor_K_CH)):
        for i in range(0, len(arrayFor_K_CH)):
            variantsPerTime[t][p] += Result[i] * eigenvectors[p][i] * math.e**(t/5 * eigenvalues[i])


for p in range(0,len(arrayFor_K_CH)):
    yForGraph = []
    for t in range(0, time):
        yForGraph.append(variantsPerTime[t][p])
    plt.plot(timeForPlot,yForGraph)
plt.xlabel('Time', color='gray')
plt.ylabel('P(t)',color='gray')
plt.legend(['P0','P1','P2','P3','P4','P5','P6','P7','P8'])
plt.show()



averageMaintenance = []
averageQueue = []
averageSystem = []
averageWait = []
averageWaitQueue = []
averageSystemTime = []
relativeThroughput = []
averageServiceTime = []
chanceBuzy = []
absoluteThroughput = []
for t in range (0,time):
    chanceOfFail = 0
    chanceBuzy.append(1 - round(variantsPerTime[t][0],5))
    relativeThroughput.append(1)
    absoluteThroughput.append(Lyambda_2 * chanceBuzy[t])

    averageMaintenance.append(0)
    for k in range(0,L):
        averageMaintenance[t] += k * round(variantsPerTime[t][k],5)
    for k in range(L,size):
        averageMaintenance[t] += L * round(variantsPerTime[t][k],5)

    averageSystem.append(0)
    for k in range(0, size):
        averageSystem[t] += k * round(variantsPerTime[t][k], 5)

    averageQueue.append(averageSystem[t] - averageMaintenance[t])

    Delta = (K - averageSystem[t]) * Lyambda  # Дельта
    averageWait.append(0)
    for k in range(0,L):
        averageWait[t] += (L-k) * round(variantsPerTime[t][k],5)

    averageWaitQueue.append((1 / Delta) * averageQueue[t])

    averageServiceTime.append((1 / Delta) * averageMaintenance[t])
    averageSystemTime.append(averageServiceTime[t] + averageWaitQueue[t])

#print(" Strain = ",systemStrain)
#print(" Strain per channel = ",systemStrainPerChannel)
#print(" Chance fail = ", chanceOfFail)
#print(" Chance Service = ", chanceOfService)
#print(" Relative throughput = ", relativeThroughput)
#print(" Absolute throughput = ", absoluteThroughput)
print(" Average maintenance = ", averageMaintenance)
plt.plot(timeForPlot,averageMaintenance)
plt.xlabel('Time', color='gray')
plt.ylabel('Nobs',color='gray')
plt.show()
print(" Average queue = ", averageQueue)
plt.plot(timeForPlot,averageQueue)
plt.xlabel('Time', color='gray')
plt.ylabel('Noch',color='gray')
plt.show()
print(" Average system = ", averageSystem)
plt.plot(timeForPlot,averageSystem)
plt.xlabel('Time', color='gray')
plt.ylabel('Nsys',color='gray')
plt.show()
print(" Average wait = ", averageWait)
plt.plot(timeForPlot,averageWait)
plt.xlabel('Time', color='gray')
plt.ylabel('Nwait',color='gray')
plt.show()
print(" Average wait Queue = ", averageWaitQueue)
plt.plot(timeForPlot,averageWaitQueue)
plt.xlabel('Time', color='gray')
plt.ylabel('Toch',color='gray')
plt.show()
print(" Average Service Time = ", averageServiceTime)
plt.plot(timeForPlot,averageServiceTime)
plt.xlabel('Time', color='gray')
plt.ylabel('Tservice',color='gray')
plt.show()
print(" Average System Time = ", averageSystemTime)
plt.plot(timeForPlot,averageSystemTime)
plt.xlabel('Time', color='gray')
plt.ylabel('Tsys',color='gray')
plt.show()