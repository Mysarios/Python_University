import math
import sympy as sp
import numpy as np
from numpy import linalg as LA

# k = 2
# lambda_ = 3.5
# mu_ = 6

k = 3
lambda_ = 2.5
mu_ = 5

# str_ = input()
# arr = []
# while str_:
#     arr.append([])
#     for el in str_.split("\t"):
#         arr[-1].append(float(el.replace(',', '.')))
#     str_ = input()

# print(arr)
arr = [[-2.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,   0.0],
       [15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -17.5, 0.0, 0.0, 2.5],
       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -15.0, 0.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -15.0, 0.0],
       [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 15.0, -15.0],
       ]

print("Система уравнений:")
p = []
uravn = []
summ_p = 0
arr_tr = np.array(arr).transpose()
for i in range(len(arr_tr)):
    p.append(sp.Symbol("p_" + str(i)))
    summ_p += p[-1]
    uravn.append(0)
    for j in range(len(arr_tr[i])):
        uravn[-1] += arr_tr[i][j] * sp.Symbol("p_" + str(j))
    print(uravn[-1], "= 0")
uravn.pop()
summ_p -= 1
uravn.append(summ_p)
p_uravn = sp.solve(uravn, p)
# print(p_uravn)


print("\nСистема Колмогорова-Чепмена для вероятностей состояний СМО в динамическом режиме")
for i in range(len(arr)):
    print("P'_" + str(i) + "(t) = ", end="")
    ind_znak = 0
    for j in range(len(arr[i])):
        if arr[i][j]:
            if ind_znak:
                print(" + ", end="")
            print(str(arr[i][j]) + "*P_" + str(j) + "(t)", end="")
            ind_znak = 1
    print()

# Этапы решения системы ОДУ Колмогорова-Чепмена
Matrix = np.array(arr)
Matrix = Matrix.transpose()
print("\nЭтапы решения системы ОДУ Колмогорова-Чепмена:")
print("Находим собственные числа и вектора матрицы")
eigenvalues, eigenvectors = LA.eig(Matrix)
for i in range(len(eigenvalues)):
    eigenvalues[i] = eigenvalues[i].round(5)
print("Собственные числа(gamma):\n", eigenvalues)
print("Собственные вектора(V):\n", eigenvectors)
print("Через систему уравнений P_i(t) = summ_j(C_j*V_ij*exp(gamma_j*t)), где i и j изменяются от 0 до k+l+1\n"
      "и условие P(0) = (1,0,...,0)T нахожим коэффициенты C_j")
c_sym_vector = []
for i in range(len(eigenvectors)):
    c_sym_vector.append(sp.Symbol("c_" + str(i)))
p_matr = []
p_matr_t_null = []
p_stach = []
for i in range(len(eigenvectors)):
    p_matr.append(0)
    for j in range(len(eigenvectors[i])):
        p_matr[-1] += c_sym_vector[j] * eigenvectors[i][j].round(3) * sp.exp(eigenvalues[j].round(3) * sp.Symbol('t'))
        if eigenvalues[j] == 0:
            p_stach.append(c_sym_vector[j] * eigenvectors[i][j].round(3))
    if not i:
        p_matr_t_null.append((p_matr[-1] - 1).subs(sp.Symbol('t'), 0))
    else:
        p_matr_t_null.append(p_matr[-1].subs(sp.Symbol('t'), 0))
c_vector = sp.solve(p_matr_t_null, c_sym_vector)
for el in c_vector:
    c_vector[el] = c_vector[el].round(5)
print(c_vector)

# Вероятности состояний СМО в динамическом режиме в аналитическом виде
print("\nВероятности состояний СМО в динамическом режиме в аналитическом виде")
i = 0
for el_p in p_matr:
    for c in c_vector:
        el_p = el_p.subs(c, c_vector[c].round(3))
    print( el_p)
    i += 1

print("\nВероятности состояний СМО в стационарном режиме в аналитическом виде")
print("Из динамического режима  Из системы уравнений")
i = 0
Sym_1 =0
Sym_2 =0
for i in range(len(p_stach)):
    for c in c_vector:
        p_stach[i] = p_stach[i].subs(c, c_vector[c])
    print(" Limit ",i," = ",p_stach[i].round(10),"  and from system = ",p_uravn[p[i]].round(10))
    Sym_1+=p_stach[i].round(10)
    Sym_2+=p_uravn[p[i]].round(10)
    i += 1
print("Check sym = ",Sym_1,"and from system =",Sym_2)
print("\nВыражения для характеристик эффективности СМО в стационарном режиме")
print("Вероятность отказа")
charakt = p_stach[-1]
print("P_отк(t) = P_(k+l) =", charakt)
print("Вероятность обслуживания")
print("P_обсл(t) = 1 - P_отк(t) =", 1 - charakt)
print("Относительная пропускная способность Q(t)")
print("Q(t) = P_обсл(t) =", 1 - charakt)
print("Абсолютная пропускная способность A(t)")
print("A(t) = lambda*Q(t) =", lambda_ * (1 - charakt))

t_ob = k/mu_
sigm_ = k**(1/2)/mu_
v = sigm_/t_ob
ro_ = lambda_/mu_
print("Среднее число заявок под обслуживанием:")
print(lambda_ * (1 - charakt) / mu_)
print("Среднее число заявок в очереди:")
# r = 0
# for i in range(k, len(p_stach)):
#     r += (i - k + 1) * p_stach[i]
# print(r)
print(ro_**2*(1+v**2)/(2*(1-ro_)))
print("Среднее время в очереди")
# t = 0
# for i in range(1, len(p_stach)):
#     t += p_stach[i]/(k*mu_)
# print(t)
print(ro_**2*(1+v**2)/(2*lambda_*(1-ro_)))
print("Среднее время обслуживания")
# t = 0
# for i in range(1, len(p_stach)):
#     t += p_stach[i]/(k*mu_)
# print(t)
print(1/mu_)
