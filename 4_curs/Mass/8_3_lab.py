import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
from numpy import linalg as LA
from sympy import *

K = 2
L = 5
size = K + L + 1
u = 0.16
l = 2.12

inputApplicationFlow = [3, 2, 3, 1, 4, 2, 4, 2, 0, 2, 1, 3, 6, 4, 2, 1, 2,
                        1, 1, 2, 2, 2, 4, 5, 1, 0, 1, 2, 0, 2, 1, 1, 3,
                        1, 1, 4, 2, 2, 3, 3, 2, 0, 4, 2, 2, 2, 1, 3, 1, 5,
                        5, 2, 2, 2, 0, 2, 1, 2, 1, 1, 3, 5, 0, 2, 2, 5, 2, 0,
                        4, 5, 1, 1, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 2, 4, 3, 2, 0, 1, 2, 2, 2, 1, 6, 3, 2, 1, 3, 1, 2, 1]
inputTimeFlow = [0.073,0.033,0.011,0.096,0.053,0.055,0.167,0.184,0.016,0.005,0.338,0.285,0.177,0.119,0.080,0.150,0.099
                ,0.158,0.262,0.589,0.008,0.124,0.239,0.312,0.085,0.362,0.134,0.035,0.134,0.333,0.152,0.018,0.019,0.136
                ,0.100,0.370,0.221,0.176,0.026,0.053,0.415,0.022,0.292,0.367,0.012,0.115,0.065,0.019,0.018,0.014,0.209
                ,0.169,0.042,0.229,0.067,0.130,0.021,0.227,0.196,0.101,0.014,0.046,0.044,0.050,0.199,0.119,0.359,0.159
                ,0.050,0.108,0.315,0.366,0.151,0.083,0.034,0.095,0.131,0.049,0.310,0.079,0.552,0.030,0.170,0.052,0.305
                ,0.181,0.466,0.038,0.092,0.316,0.020,0.096,0.320,0.171,0.067,0.109,0.869,0.207,0.433,0.027]
def replace_string(s):
    buf =""
    for c in s:
        if c == ",":
            buf +="."
            continue
        if c == " ":
            buf += ","
            continue
        buf += c
    print(buf)

Lyambda = np.mean(inputApplicationFlow)
print(Lyambda)


vals = expon.cdf(inputTimeFlow)
print("Принимает гипотезу = ", np.allclose(inputTimeFlow, expon.ppf(vals)))

vals = poisson.cdf(inputApplicationFlow, Lyambda)
print("Принимает гипотезу = ",np.allclose(inputApplicationFlow, poisson.ppf(vals, Lyambda)))

Lyambda_2 = np.mean(inputTimeFlow)
print(Lyambda_2)

Array_intens = [0]*(size)
for i in range (0,len(Array_intens)):
    Array_intens[i] = [0]*(size)


#Uu = Symbol('u')
intens = 1
for low in range(1,len(Array_intens)):
    Array_intens[low][low-1] = intens * u
    if intens < K:
        intens +=1

#Ll = Symbol('l')
for high in range(0,len(Array_intens)-1):
    Array_intens[high][high+1] = l


for middle in range(0,len(Array_intens)-1):
    Array_intens[middle][middle] -= (Array_intens[middle][middle-1] +Array_intens[middle][middle+1])

Array_intens[len(Array_intens)-1][len(Array_intens)-1] = -K*u

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


for i in range(0, len(Static_eval)):
    print(Static_eval[i])
#for i in range(0, len(Static_eval)):
    #print(arrayFor_K_CH[i])

Result = 0
for solution in linsolve(Static_eval, C_symbols):
    Result = solution

print("C coefs  =",Result)


Transp_Array_intens = Array_intens.copy()
Transp_Array_intens = np.array(Transp_Array_intens)
Transp_Array_intens = Transp_Array_intens.transpose()
for i in range(0,K+L):
    print(Transp_Array_intens[i])

print("System   = ")
for i in range(0, size):
    buf = 0
    for j in range(0, size):
        arrayFor_K_CH[i] += Result[j] * eigenvectors[i][j] * math.e**(Tt * eigenvalues[j])
        buf += round(Result[j] * eigenvectors[i][j],5) * round(math.e,4)**(Tt * round(eigenvalues[j],4))
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

averageMaintenance = 0
for k in range(0,K):
    averageMaintenance += k * round(limit_expr[k],5)
for k in range(K,size):
    averageMaintenance += K * round(limit_expr[k],5)

print(" Average maintenance = ", averageMaintenance)

averageQueue = 0
for k in range(K+1,size):
    averageQueue += (k-K) * round(limit_expr[k],5)
print(" Average queue = ", averageQueue)

averageSystem = averageQueue + averageMaintenance
print(" Average system = ", averageSystem)

averageWait = 0
for k in range(0,K):
    averageWait += (K-k) * round(limit_expr[k],5)
print(" Average wait = ", averageWait)

averageWaitQueue_Time = 0
j = 1
for k in range(K,size-1):
    print(j," and k=",k)
    averageWaitQueue_Time += (j/(K*Lyambda_2)) * round(limit_expr[k],5)
    j+=1

print(" Average wait Queue Time = ", averageWaitQueue_Time) #&&

averageServiceTime = relativeThroughput/Lyambda_2
print(" Average Service Time = ", averageServiceTime)

averageSystemTime = averageServiceTime + averageWaitQueue_Time
print(" Average System Time = ", averageSystemTime)

print("Formuls Littl:")

T_Servise = averageQueue/Lyambda #&&
print(" T_wait = ", T_Servise) #&&

T_Maintenance = averageMaintenance / Lyambda
print(" T_Service Time = ", T_Maintenance)
T_System = averageSystem/Lyambda
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