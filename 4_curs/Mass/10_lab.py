import math

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import *
from numpy import linalg as LA
from sympy import *

Queue = 2
Service = 1
size = Queue + Service + 1
u = 6.25
l = 2.12
w = l * 0.1

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
Lyambda_2 = 1/np.mean(inputTimeFlow)

vals = expon.cdf(inputTimeFlow)
print("Принимает гипотезу о соответствии первого распределения распределению Пуасона= ", np.allclose(inputTimeFlow, expon.ppf(vals)))
print("Интенсивность входного потока = ",Lyambda)

vals = poisson.cdf(inputApplicationFlow, Lyambda)
print("Принимает гипотезу о соответствии первого распределения распределению экспонентациальному = ",np.allclose(inputApplicationFlow, poisson.ppf(vals, Lyambda)))
print("Интенсивность выходного потока = ",Lyambda_2)

u = Lyambda_2
l = Lyambda

#Получение первого для очереди = 5
print("Queue = 5")
#Создаем матрицу А
Array_intens = [0]*(size)
for i in range (0,len(Array_intens)):
    Array_intens[i] = [0]*(size)
intens = 0
#Нижняя линия
for low in range(1,len(Array_intens)):
    Array_intens[low][low-1] = u + intens * w
    intens +=1
#Верхняя линия
for high in range(0,len(Array_intens)-1):
    Array_intens[high][high+1] = l
#Средняя линия
for middle in range(0,len(Array_intens)-1):
    Array_intens[middle][middle] -= (Array_intens[middle][middle-1] +Array_intens[middle][middle+1])

Array_intens[len(Array_intens)-1][len(Array_intens)-1] = -Queue * w - u

print("base =")
for i in range(0,size):
    print(Array_intens[i])

#Транспонируем матрицу А
Matrix = np.array(Array_intens)
Matrix = Matrix.transpose()
print("transp base = =")
for i in range(0,size):
    print(Matrix[i])

#Ищем собственные вектора и числа
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

#Создание системы для К-Ч
C_symbols =[0]*(size)
P_symbols = [0]*(size)
arrayFor_K_CH = [0]*(size)
Static_eval = [0]*(size)
buf = [0]*(size)
for i in range (0,size):
    C_symbols[i] = Symbol('C' + str(i))
    P_symbols[i] = Symbol('P' + str(i))
arrayFor_K_CH[0] = 0
Static_eval[0] = -1
buf[0] = -1

Tt = Symbol('t')
for i in range(0, size):
    for j in range(0,size):
        Static_eval[i] += C_symbols[j] * eigenvectors[i][j]
        buf[i] += C_symbols[j] * round(eigenvectors[i][j],4)
for i in range(0, len(Static_eval)):
    print(buf[i])

print("System with C  = ")
for i in range(0, size):
    buf = 0
    for j in range(0, size):
        buf += C_symbols[j] * round(eigenvectors[i][j],2) * Symbol('e')**(Tt * round(eigenvalues[j],1))
    print(buf)

#Поиск C коэфы
Result = 0
for solution in linsolve(Static_eval, C_symbols):
    Result = solution
for c in Result:
    print("C coefs  =",round(c,3))


Transp_Array_intens = Array_intens.copy()
Transp_Array_intens = np.array(Transp_Array_intens)
Transp_Array_intens = Transp_Array_intens.transpose()
for i in range(0, Queue + Service):
    print(Transp_Array_intens[i])

#Вывод системы
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

last_result = limit_expr
#for m in range(3,20):
check_eps = 1
while check_eps > 0.00001:
    Queue += 1
    size = Queue + Service + 1
    #Создаем матрицу А
    Array_intens = [0]*(size)
    for i in range (0,len(Array_intens)):
        Array_intens[i] = [0]*(size)
    intens = 1
    # Нижняя линия
    for low in range(1, len(Array_intens)):
        Array_intens[low][low - 1] = u + intens * w
        intens += 1
    # Верхняя линия
    for high in range(0, len(Array_intens) - 1):
        Array_intens[high][high + 1] = l
    # Средняя линия
    for middle in range(0, len(Array_intens) - 1):
        Array_intens[middle][middle] -= (Array_intens[middle][middle - 1] + Array_intens[middle][middle + 1])
    Array_intens[len(Array_intens) - 1][len(Array_intens) - 1] = -Queue * w - u
    #print("base =")
    #for i in range(0,size):
        #print(Array_intens[i])

    #Транспонируем матрицу А
    Matrix = np.array(Array_intens)
    Matrix = Matrix.transpose()
    #print("transp base = =")
    #for i in range(0,size):
        #print(Matrix[i])

    #Ищем собственные вектора и числа
    eigenvalues, eigenvectors = LA.eig(Matrix)
    #print("E values = ")
    for i in range(0,size):
        if abs(eigenvalues[i]) - 0.005 < 0:
            eigenvalues[i] = 0
        #print("E" + str(i) + " ",eigenvalues[i])
    #print("E vectors = ")
    for i in range(0, size):
        buf = "EV"+ str(i)
        for j in range(0, size):
            buf += ": {:.4f}".format(eigenvectors[i][j]) + " "
        #print(buf)

    #Создание системы для К-Ч
    C_symbols =[0]*(size)
    P_symbols = [0]*(size)
    arrayFor_K_CH = [0]*(size)
    Static_eval = [0]*(size)
    buf = [0]*(size)
    for i in range (0,size):
        C_symbols[i] = Symbol('C' + str(i))
        P_symbols[i] = Symbol('P' + str(i))
    arrayFor_K_CH[0] = 0
    Static_eval[0] = -1
    buf[0] = -1

    Tt = Symbol('t')
    for i in range(0, size):
        for j in range(0,size):
            Static_eval[i] += C_symbols[j] * eigenvectors[i][j]
            buf[i] += C_symbols[j] * round(eigenvectors[i][j],4)
    #for i in range(0, len(Static_eval)):
        #print(buf[i])

    #Поиск C коэфы
    Result = 0
    for solution in linsolve(Static_eval, C_symbols):
        Result = solution
    #for c in Result:
        #print("C coefs  =",round(c,3))


    Transp_Array_intens = Array_intens.copy()
    Transp_Array_intens = np.array(Transp_Array_intens)
    Transp_Array_intens = Transp_Array_intens.transpose()
    #for i in range(0, Queue + Service):
        #print(Transp_Array_intens[i])

    #Вывод системы
    #print("System   = ")
    for i in range(0, size):
        buf = 0
        for j in range(0, size):
            arrayFor_K_CH[i] += Result[j] * eigenvectors[i][j] * math.e**(Tt * eigenvalues[j])
            buf += round(Result[j] * eigenvectors[i][j],5) * round(math.e,2)**(Tt * round(eigenvalues[j],2))
        #print(buf)
    print("Results(",Queue,")   = ")
    Sym = 0
    limit_expr = []
    for i in range (0,len(arrayFor_K_CH)):
        limit_expr.append(limit(arrayFor_K_CH[i], Tt, oo))
        Sym += limit_expr[i]
        print(" Limit ",i," = ",round(limit_expr[i],9))
    print("Check sym =", round(Sym,2))
    check_eps = 0
    for result_index in range(0,len(last_result)):
        check_eps+= round(abs(last_result[result_index] - limit_expr[result_index]),5)
    print(check_eps)
    last_result = limit_expr

#Вывод характеристик
print("Effectivnes params by t = oo:")
systemStrain = Lyambda/Lyambda_2
print(" Strain = ",systemStrain)
systemStrainPerChannel = systemStrain/2
print(" Strain per channel = ",systemStrainPerChannel)
chanceOfFail = 0
print(" Chance fail = ", chanceOfFail)
chanceOfService = 1
print(" Chance Service = ", chanceOfService)

averageQueue = 0                                                                          #Среднее заявок  в очереди Nоч
for k in range(1, size):
    averageQueue += k * round(limit_expr[k], 5)
print(" Average queue = ", averageQueue)

absoluteThroughput = l - w * averageQueue
print(" Absolute throughput = ", absoluteThroughput)
relativeThroughput = absoluteThroughput/u
print(" Relative throughput = ", relativeThroughput)


averageMaintenance = absoluteThroughput/ u                                        #Среднее заявок под обслуживанием  Nоб
print(" Average maintenance = ", averageMaintenance)

averageSystem = averageQueue + averageMaintenance                                       #Среднее заявок в системе  Nсист
print(" Average system = ", averageSystem)

averageWait = 0                                                                                     #Среднее простоя Nпр
for k in range(0, Queue):
    averageWait += (Queue - k) * round(limit_expr[k], 5)
print(" Average wait = ", averageWait)

averageWaitQueue_Time = 0                                                                   #Среднее время в очереди Tоч
for k in range(1, size):
    averageWaitQueue_Time += (1/w) * round(limit_expr[k], 5)
print(" Average wait Queue Time = ", averageWaitQueue_Time)                                                          #&&

averageServiceTime = l                                                                   #Среднее время обслуживания Tоб
print(" Average Service Time = ", averageServiceTime)

averageSystemTime = averageServiceTime + averageWaitQueue_Time                            #Среднее время в системе Tсист
print(" Average System Time = ", averageSystemTime)



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
absoluteThroughput = []
averageServiceTime = []
for t in range (0,time):
    chanceOfFail = 0
    chanceOfService = 1

    averageQueue.append(0)
    for k in range(1, size):
        averageQueue[t] += k * round(variantsPerTime[t][k], 5)

    absoluteThroughput.append(l - w * averageQueue[t])
    relativeThroughput.append(absoluteThroughput[t]/u)
    averageMaintenance.append(absoluteThroughput[t]/u)

    averageSystem.append(averageQueue[t] + averageMaintenance[t])

    averageWait.append(0)
    for k in range(0, Queue):
        averageWait[t] += (Queue - k) * round(variantsPerTime[t][k], 5)

    averageWaitQueue.append(0)
    for k in range(1, size):
        averageWaitQueue[t] += (1/w) * round(variantsPerTime[t][k], 5)

    averageServiceTime.append(l)
    averageSystemTime.append(averageServiceTime[t] + averageWaitQueue[t])


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