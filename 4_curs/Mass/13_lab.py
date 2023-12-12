import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
from numpy import linalg as LA
from sympy import *

r = 3
L = 10
size = r + L + 1
u = 5
l = 2.5

Lyambda = l
Lyambda_2 = u

u = Lyambda_2
l = Lyambda

Array_intens = [0]*(size)
for i in range (0,len(Array_intens)):
    Array_intens[i] = [0]*(size)


#Uu = Symbol('u')
intens = 1
for low in range(1,len(Array_intens)):
    Array_intens[low][low-1] = r*u

#Ll = Symbol('l')
for high in range(0,len(Array_intens)-r):
    Array_intens[high][high+r] = l


for middle in range(0,len(Array_intens)-r):
    Array_intens[middle][middle] -= (Array_intens[middle][middle-1] +Array_intens[middle][middle+r])
for i in range(1,r+1):
    Array_intens[len(Array_intens)-i][len(Array_intens)-i] = -r*u

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
    print("E" + str(i) + " {:.3f}".format(eigenvalues[i]))

print("E vectors = ")

for i in range(0, size):
    buf = "EV"+ str(i)
    for j in range(0, size):
        buf += ": {:.3f}".format(eigenvectors[i][j]) + " "
    print(buf)


C_symbols =[0]*(size)
P_symbols = [0]*(size)
arrayFor_K_CH = [0]*(size)
Static_eval = [0]*(size)
buf = [0]*(size)
for i in range (0,size):
    #arrayFor_K_CH[i] = [0]*9
    #Static_eval[i] = [0]*9
    C_symbols[i] = Symbol('C' + str(i))
    P_symbols[i] = Symbol('P' + str(i))

arrayFor_K_CH[0] = 0
Static_eval[0] = -1
buf[0] = -1

Tt = Symbol('t')
for i in range(0, size):
    for j in range(0,size):
        #arrayFor_K_CH[i] += C_symbols[j] * eigenvectors[i][j] * math.e**(Tt * eigenvalues[j])
        Static_eval[i] += C_symbols[j] * eigenvectors[i][j]
        buf[i] += C_symbols[j] * round(eigenvectors[i][j],4)


for i in range(0, len(Static_eval)):
    print(buf[i])
#for i in range(0, len(Static_eval)):
    #print(arrayFor_K_CH[i])

Result = 0
for solution in linsolve(Static_eval, C_symbols):
    Result = solution
for c in Result:
    print("C coefs  =",round(c,3))


Transp_Array_intens = Array_intens.copy()
Transp_Array_intens = np.array(Transp_Array_intens)
Transp_Array_intens = Transp_Array_intens.transpose()

print("System   = ")
for i in range(0, size):
    buf = 0
    for j in range(0, size):
        arrayFor_K_CH[i] += Result[j] * eigenvectors[i][j] * math.e**(Tt * eigenvalues[j])
        buf += Result[j] * eigenvectors[i][j] * round(math.e,2)**(Tt * eigenvalues[j],2)
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

chanceOfFail = round(limit_expr[len(limit_expr)-1],5)
print(" Chance fail = ", chanceOfFail)
chanceOfService = 1 - chanceOfFail
print(" Chance Service = ", chanceOfService)
relativeThroughput = chanceOfService
print(" Relative throughput = ", relativeThroughput)
absoluteThroughput = Lyambda * chanceOfService
print(" Absolute throughput = ", absoluteThroughput)

averageMaintenance = 0                                  #Среднее заявок под обслуживанием  Nоб
for k in range(0,K):
    averageMaintenance += k * round(limit_expr[k],5)
for k in range(K,size):
    averageMaintenance += K * round(limit_expr[k],5)

print(" Average maintenance = ", averageMaintenance)

averageQueue = 0                                        #Среднее заявок  в очереди Nоч
for k in range(K+1,size):
    averageQueue += (k-K) * round(limit_expr[k],5)
print(" Average queue = ", averageQueue)

averageSystem = averageQueue + averageMaintenance       #Среднее пасивов  Nсист
print(" Average system = ", averageSystem)

averageWait = 0                                         #Среднее простоя Nпр
for k in range(0,K):
    averageWait += (K-k) * round(limit_expr[k],5)
print(" Average wait = ", averageWait)

averageWaitQueue_Time = 0                               #Среднее время в очереди Tоч
j = 1
for k in range(K,size-1):
    print(j," and k=",k)
    averageWaitQueue_Time += (j/(K*Lyambda_2)) * round(limit_expr[k],5)
    j+=1

print(" Average wait Queue Time = ", averageWaitQueue_Time) #&&

averageServiceTime = relativeThroughput/Lyambda_2       #Среднее время обслуживания Tоб
print(" Average Service Time = ", averageServiceTime)

averageSystemTime = averageServiceTime + averageWaitQueue_Time  #Среднее время в системе Tсист
print(" Average System Time = ", averageSystemTime)

print("Formuls Littl:")

T_Maintenance = averageQueue/Lyambda #&&                        #Среднее время в очереди Tоч
print(" T_wait = ", T_Maintenance) #&&

T_Servise = averageMaintenance / Lyambda                #Среднее время обслуживания Tоб
print(" T_Service Time = ", T_Servise)
T_System = averageSystem/Lyambda                            #Среднее время в системе Tсист
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
plt.legend(['P0','P1','P2','P3','P4','P5','P6','P7'])
plt.show()



averageMaintenance = []
averageQueue = []
averageSystem = []
averageWait = []
averageWaitQueue = []
averageSystemTime = []
relativeThroughput = []
averageServiceTime = []
for t in range (0,time):
    chanceOfFail = round(variantsPerTime[t][len(limit_expr) - 1], 5)
    chanceOfService = 1 - chanceOfFail
    relativeThroughput.append(chanceOfService)

    averageMaintenance.append(0)
    for k in range(0,K):
        averageMaintenance[t] += k * round(variantsPerTime[t][k],5)
    for k in range(K,size):
        averageMaintenance[t] += K * round(variantsPerTime[t][k],5)

    averageQueue.append(0)
    for k in range(K+1,size):
        averageQueue[t] += (k-K) * round(variantsPerTime[t][k],5)

    averageSystem.append(averageQueue[t] + averageMaintenance[t])

    averageWait.append(0)
    for k in range(0,K):
        averageWait[t] += (K-k) * round(variantsPerTime[t][k],5)

    averageWaitQueue.append(0)
    j = 1
    for k in range(K,size-1):
        averageWaitQueue[t] += (j/(K*Lyambda_2)) * round(variantsPerTime[t][k],5)
        j += 1

    averageServiceTime.append(relativeThroughput[t]/Lyambda_2)
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