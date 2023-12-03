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

inputType_1 = [138.1,90.5,43.46,96.36,42.29,131.51,13.58,22.33,29.73,27.61,25.61,40.66,102.74,63.02,91.12,30.69,10,15.4
                ,24.66,30.6,123.28,72.26,86.63,23.69,57.85,87.84,88.61,119.68,121.76,65.46,60.51,94.79,79.19,57.1,56.06
                ,66.31,86.75,98.59,81.08,99.85,65.44,92.97,42.63,15.25,21.83,26.53,60.93,30.28,80.91,45.24,114.13,13.41
                ,108.39,59.15,41.43,22.32,15.99,110.09,71.61,46.88,48.91,29.44,126.09,79.41,93.91,77.92,87.92,31.22
                ,153.47,97.04,17.09,111.93,48.32,151.48,90.93,22.72,39.05,45.38,96.19,72.36,73.91,31.53,21.08,146.47,44.21
                ,42.83,77.39,19.48,29.84,71.97,34.95,7.64,29.05,91.93,83.23,48.1,58.79,37.17,55.49,158.7,121.15,324.83,313.72,348.98,334.01
                ,349.55,292,299.18,348.94,320.83,252.23,248.28,378.61,284.08,333.73,353.32,309.83,379.84,252.85,377.14,265.25,355.06,334.92
                ,374.16,229.23,298.67,237.3,362.43,297.43,319.55,284.54,311.17,308.87,314.33,273.92,296.71,319.66,319.44,210.03,335.32,247.46
                ,292.16,273.16,313.6,305.69,307.94,283.84,333.78,360.36,287.33,252.3,346.39,322.13,330.87,264.32,241.42,305.45,262.8,253.93
                ,259.75,349.07,326.61,280.34,352.4,312.94,291.8,304.08,309.71,305.24,224.25,309.86,273.45,348,273.39,292.28,197.68,328.11
                ,253.8,302.07,339.75,330.78,250.85,271.81,333.14,310.55,252.12,325.83,340.24,276.18,318.41,286.07,307.14,309.19,288.52,348.5,280.88]
inputType_2 = [4.71,2.58,2.82,2.69,3.02,4.19,3.02,3.79,4.26,4.51,4.12,5.2,4.76,2.87,3.55,3.86,2.06,4.88,3.58,2.94,5.63
                ,3.16,6.1,6.34,3.99,2.62,3.78,3.79,2.82,3.53,6.28,2.82,5.59,3.68,3.19,3.57,4.15,2.92,5.79,3.27,3.76,5.08,3.57,4.37
                ,3.53,4.3,4.7,4.39,3.39,3.42,5,3.64,4.81,3.48,4.14,3.57,3.03,3.45,3,2.04,5.25,4.72,4.23,4.51,3.89,5.72,4.29,2.75,5.34
                ,2.97,4.02,4.17,3.68,3.65,4.23,4.03,5.3,5.09,3.54,3.44,4.77,3.64,3.53,3.32,4.83,3.6,1.83,5.29,4.05,3.77,3.42,3.04,2.49
                ,2,4.74,3.98,3.89,4.26,3.36,5.22,2.35,3.62,4.88,4.37,2.64,4.63,6.03,4.8,5.27,6.22,5.17,1.82,4.8,4.21,3.75,3.06,2.87,2.64
                ,3.55,2.39,4.54,3.29,4.35,4.08,4.95,4.65,4.98,4.94,3.36,4.29,3.94,5.49,3.92,4.06,4.78,5.9,3.37,2.26,4.14,4.11,3.57,2.99
                ,3.14,4.8,4.5,2.09,2.85,3.42,4.27,3.55,3.62,3.85,2.19,2.76,3.97,4.2,2.32,4.25,1.57,3.44,3.58,4.95,1.48,3.39,3.51,3.94,3.31
                ,3.58,4.83,5.15,4.78,4.34,5.27,4.16,5.5,5.01,2.19,3.46,5.49,3.99,4.54,5.22,5.11,6.42,5.85,5.15,4,4.14,3.26,3.48,3.77,4.68
                ,3.52,4.3,6.32,4.6]
inputType_3 = [12.57,2.77,3.75,0.56,6.05,0.07,0.59,2.15,0.78,8.61,0.59,6.1,5.76,0.23,2.41,1.53,6.52,1.04,6.27,3.47,3.41
                ,2.1,1.85,1.65,0.67,0.6,1.53,0.88,3.66,2.44,13.11,16.72,7.51,3.11,6.1,0.44,7.07,5.26,7.14,7.15,0.24,1.26,11.19,0.88
                ,4.57,0.39,1.06,7.9,9.13,2.54,2.19,0.8,2.73,15,3.96,15,1.08,2.72,6.13,3.73,1.53,9.49,0.78,3.84,2.7,14.02,8.55,3.28
                ,6.74,0.26,5.09,1.25,4.97,10.12,1.73,0.7,0.04,15,4.86,0.08,1.8,11.58,3.93,0.91,3.36,0.86,2.97,5.13,2.34,1.95,0.12
                ,6.08,0.17,2.61,5.54,0.47,10.44,4.85,6.36,6.93,15,9.15,3.14,1.54,0.91,1.04,1.55,1.68,8.32,8.01,3.88,7.5,2.76,0.67
                ,1.9,3.6,2.77,2.67,15,8.1,6,9.46,8.11,5.5,4.76,3.96,2.9,15.78,9.15,0.79,5.23,1.63,5.54,3.3,4.7,1.45,1.29,7.2,8.49
                ,5.06,1.95,2.43,4.17,4,7.97,2.46,14.73,2.81,2.6,11.61,10.1,2.39,2.66,4.01,2.37,4.33,15,1.93,8.22,8.79,6.27,3.07
                ,6.28,5.77,15,2.61,5.27,3,2.84,8.05,3.78,12.15,15,0.45,15,3.95,8.26,1.72,4.27,2.87,14.4,0.82,4.34,10.84,3.88
                ,11.78,1.51,5.82,18.61,8.7,1.58,3.65,8.47,5.58,16.57,4.44]

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
Lyambda_2 = 1/np.mean(inputTimeFlow)

vals = expon.cdf(inputTimeFlow)
print("Принимает гипотезу о соответствии первого распределения распределению Пуасона= ", np.allclose(inputTimeFlow, expon.ppf(vals)))
print("Интенсивность входного потока = ",Lyambda)

vals = poisson.cdf(inputApplicationFlow, Lyambda)
print("Принимает гипотезу о соответствии первого распределения распределению экспонентациальному = ",np.allclose(inputApplicationFlow, poisson.ppf(vals, Lyambda)))
print("Интенсивность выходного потока = ",Lyambda_2)

u = Lyambda_2
l = Lyambda

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
for i in range(0,K+L):
    print(Transp_Array_intens[i])

print("System   = ")
for i in range(0, size):
    buf = 0
    for j in range(0, size):
        arrayFor_K_CH[i] += Result[j] * eigenvectors[i][j] * math.e**(Tt * eigenvalues[j])
        buf += round(Result[j] * eigenvectors[i][j],5) * round(math.e,2)**(Tt * round(eigenvalues[j],2))
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