# 1)	  - распределение Коши с плотностью  ,  , где параметр   - неизвестен (оценку параметра   найти численно, путем минимизации значения статистики критерия хи-квадрат);
# 2)	  - нормальный закон распределения с плотностью   ,  , где параметры   и   - неизвестны.
import math
from matplotlib import pyplot as plt
import numpy as np
from sympy import *
import sympy as sym

Sigma = Symbol('Sigma')
A = Symbol('a')
X = Symbol('x')
t = Symbol('t')

PI = 3.141592653589793

def Get_M(Array):
    M=Array[0]
    for i in range (1,len(Array)):
        M+=Array[i]

    M/=len(Array)
    return M
def Get_D(Array):
    D = np.var(Array)
    D = D**(1/2.0)
    return D
def Get_Function_Normal_Raspredelenia(x,a,sigma):
    Buf_0 = 1.0/( sigma * (2*math.pi)**(1/2.0))
    Buf_1 = math.exp(-((x - a)**2) / (2 * (sigma**2)))
    Function = Buf_0 * Buf_1

    return Function
def Create_Symb_Function_Normal(sigma,a):
    Function = ( sym.exp(( -((t - a) ** 2) / (2 * (sigma ** 2)) ) ))/ (sigma * (2 * PI) ** (1 / 2.0))
    return Function.subs([(A,a),(Sigma,sigma)])
def Get_integral(Function,x,a):
    if x > a:
        Result = integrate(Function, (t, a, x))
    if x <= a:
        Result = integrate(Function, (t, x, a))

    Result.subs(t, x)
    Result = 1/2 + Result
    return Result.evalf()

Base_Array = [0.1, -1.53, -0.94, 0.21, 0.77, 1.10, 0.23, -0.15, 0.79, -0.71,
              1.17, 0.01, 0.45, 1.55, 1.48, -0.09, 0.01, 1.00, 1.25, 1.35,
              0.52, -1.61, 2.16, 0.64, 0.19, 0.02, 0.20, 1.43, 0.74, -0.21,
              0.41, 0.80, -0.41, 0.31, -1.26, 0.75, 1.05, 2.04, -0.42, -1.06,
              0.33, -0.30, -0.34, -0.10, -1.54, 0.67, -0.40, -0.15, 0.98, -1.04,
              1.55, -1.58, 1.78, -0.71, 0.75, 0.48, -0.18, 0.49, -0.07, 0.90,
              1.04, 2.75, 1.03, 0.76, -2.53, 0.27, 0.92, -1.17, -0.85, -1.83,
              -0.35, -1.07, -0.02, 1.64, 0.35, -0.86, -0.06, 0.69, 2.16, -0.54,
              1.20, -0.57, 1.57, -0.05, 0.34, 0.83, -0.28, 0.48, 1.85, 0.93,
              0.91, -1.50, -1.08, 0.53, -0.53, 0.29, 0.77, -1.13, -0.76, 2.30
              ]
Base_Array_Sorted = sorted(Base_Array)
n = len(Base_Array)
Min_num = min(Base_Array)
Max_num = max(Base_Array)


# Частота попадания данных в период
Count_Intervals = int(1 + math.log(n,2))
Absolut_max = abs(Min_num)+Max_num

Interval_step = Absolut_max/Count_Intervals
periodicity= []

Count = 0
Now_max = Min_num + Interval_step
for i in range(0,n):
    if Base_Array_Sorted[i] < Now_max:
        Count += 1
    else:
        periodicity.append(Count)
        Count=1
        Now_max+=Interval_step
#print(periodicity)



# Неподходит - меняем период
Count_Intervals = int(1 + math.log(n,2))
Absolut_max = abs(Min_num)+Max_num

Interval_step = Absolut_max/(Count_Intervals-2)
periodicity= []
Middle = []
relative_periodicity = []
solidity_relative_periodicity = []
Count = 0

Now_max = Min_num + Interval_step
Middle.append( (Min_num + Min_num + Interval_step) / 2)
for i in range(0,n):
    if Base_Array_Sorted[i] < Now_max:
        Count += 1
    else:
        periodicity.append(Count)
        Count=1
        Middle.append( (Now_max + Now_max + Interval_step ) / 2)
        Now_max+=Interval_step
Middle = Middle[:-1]
print("Period: " ,periodicity)

for i in range(0 , len(periodicity)):
    relative_periodicity.append( periodicity[i]/ len(Base_Array_Sorted) )
    solidity_relative_periodicity.append( relative_periodicity[i]/Interval_step )

print("Miidle num: " , Middle)
print("Realitive_per : " , relative_periodicity)
print("Solidy: " ,solidity_relative_periodicity)
# Вроде норм
# Graph :

fig, ax = plt.subplots()
ax.hist(Base_Array_Sorted,Middle)
plt.show()

#Коши:

#Нормальный:
a = Get_M(Base_Array)
sigma =Get_D(Base_Array)
print("a = ",a," Sigma =",sigma)

solidity_Teoretic_periodicity = []
for i in range (0, len(Middle)):
    solidity_Teoretic_periodicity.append(Get_Function_Normal_Raspredelenia(Middle[i],a,sigma))
print()
print("Per: ", solidity_relative_periodicity)
print("Teor per: ", solidity_Teoretic_periodicity)
print()

results_normal = []
results_normal_last = []
Normal_function = Create_Symb_Function_Normal(sigma,a)

check = []
New_number_to_Integral = Min_num
for i in range(1, len(Middle)):
    check.append(New_number_to_Integral)
    results_normal.append(Get_integral(Normal_function, New_number_to_Integral + Interval_step,a))
    results_normal_last.append(Get_integral(Normal_function, New_number_to_Integral,a))
    New_number_to_Integral +=Interval_step

results_normal.append(Get_integral(Normal_function, np.inf,a))
results_normal_last.append(Get_integral(Normal_function, New_number_to_Integral,a))
check.append(New_number_to_Integral)

print()
print(check)
print(results_normal)
print(results_normal_last)
print()

p_i = []
for i in range(0, len(results_normal)):
    p_i.append(results_normal[i] - results_normal_last[i])

print("Veroiatnosti : ", p_i)
Wait_nums = []
Check_pi = 0
for i in range(0, len(p_i)):
    Wait_nums.append(p_i[i] * n)
    Check_pi += p_i[i]
print( Check_pi)
print("Wait nums : ", Wait_nums)
#New
Buf_num = []
for i in range(0, len(p_i)):
    Buf_num.append( ((periodicity[i] - abs(Wait_nums[i]) )**2) / Wait_nums[i] )

print("per : ", periodicity)
print("New nums : ", Buf_num)

Pirson_num = 0
for i in range(0, len(p_i)):
    Pirson_num += Buf_num[i]
print("Pirson_num :",Pirson_num)
print(" Teor num :", 11.07)