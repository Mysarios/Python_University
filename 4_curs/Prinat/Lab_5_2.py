from sympy import *

Inp_type:int = int(input("Введите 0 для задания своей матрицы или 1 для тестовой(по варианту) или 2 для тестовой смешанных: "))
size_A =0
size_B =0

if Inp_type == 0:
    size_A:int = int(input("Введите кол-во стратегий для А: "))
    size_B:int = int(input("Введите кол-во стратегий для B: "))
    matrix_a =[]
    for i in range(0,size_A):
        matrix_a.append([])
        for j in range(0,size_B):
            print("Введите элемент A(" + str(i+1) + "," + str(j+1) + ")")
            matrix_a[i].append(int(input("= ")))
if Inp_type == 1:
    size_A:int = 5
    size_B:int = 5
    matrix_a =[[-3,0,0,5,1],[2,1,7,6,3],[2,-5,4,7,8],[5,0,-3,1,4],[0,1,6,7,-1]]
if Inp_type == 2:
    size_A:int = 3
    size_B:int = 3
    matrix_a =[[4,7,2],[7,3,2],[2,1,8]]

min_A = []*size_A
max_B = []*size_B

#Find down_line
for i in range(0,size_A):
    min_val = matrix_a[i][0]
    for j in range(0,size_B):
        if matrix_a[i][j] < min_val:
            min_val = matrix_a[i][j]
    min_A.append(min_val)

print("Max string A = ",min_A)
print("Нижняя граница игры = ",max(min_A))

# Find up_line
for i in range(0, size_B):
    max_val = matrix_a[0][i]
    for j in range(0, size_A):
        if matrix_a[j][i] > max_val:
            max_val = matrix_a[j][i]
    max_B.append(max_val)

print("Min string B = ",max_B)
print("Верхняя границы игры = ", min(max_B))

if min(max_B) == max(min_A) :
    index_a = 0
    for i in range(0,len(min_A)):
        if max(min_A) == min_A[i]:
            index_a = i
            break
    index_b = 0
    for i in range(0, len(max_B)):
        if min(max_B) == max_B[i]:
            index_b = i
            break
    print("Седловая точка существует в A(",index_a+1,") B(",index_b+1,")")
else:
    print("Седловая точка не существует")

print("Вывод матрицы после доминирования: ")

a_ok = 1
b_ok = 1
new_size_a = size_A
new_size_b = size_B

new_matrix = matrix_a
print("Start matrix =",new_matrix)
while a_ok or b_ok:
    print("Start Dom")
    #print("Size A =",new_size_a)
    #print("Size B =",new_size_b)

    a_ok = 0
    b_ok = 0
    index_a_for_delete = set()


    for start_row in range(0, new_size_a-1):
        for index_a in range(start_row+1,new_size_a):
            index_buf = index_a
            for index_b in range(0,new_size_b):
                if new_matrix[index_a][index_b] > new_matrix[start_row][index_b]:
                    index_buf = -1
                    break
            if index_buf != -1:
                a_ok = 1
                index_a_for_delete.add(index_buf)

    for start_row in range(new_size_a-1, 0,-1):
        for index_a in range(start_row-1,-1,-1):
            index_buf = index_a
            for index_b in range(0,new_size_b):
                if new_matrix[index_a][index_b] > new_matrix[start_row][index_b]:
                    index_buf = -1
                    break
            if index_buf != -1:
                a_ok =1
                index_a_for_delete.add(index_buf)

    new_matrix_buf = []*len(index_a_for_delete)
    new_size_a_buf = 0
    for i in range(0,new_size_a):
        if not (i in index_a_for_delete):
            new_size_a_buf += 1
            new_matrix_buf.append(new_matrix[i])
    new_size_a = new_size_a_buf
    new_matrix = new_matrix_buf
    index_b_for_delete = set()

    for start_column in range(0, new_size_b - 1):
        for index_b in range(start_column + 1, new_size_b):
            index_buf = index_b
            for index_a in range(0, new_size_a):
                if new_matrix[index_a][index_b] < new_matrix[index_a][start_column]:
                    index_buf = -1
                    break
            if index_buf != -1:
                b_ok = 1
                index_b_for_delete.add(index_buf)
    for start_column in range(new_size_b - 1, 0, -1):
        for index_b in range(start_column - 1, -1, -1):
            index_buf = index_b
            for index_a in range(0, new_size_a):
                if new_matrix[index_a][index_b] < new_matrix[index_a][start_column]:
                    index_buf = -1
                    break
            if index_buf != -1:
                b_ok = 1
                index_b_for_delete.add(index_buf)


    new_matrix_buf = []
    new_size_b_buf = 0
    for index_a in range(0,new_size_a):
        new_matrix_buf.append([])
        for index_b in range(0, new_size_b):
            if not (index_b in index_b_for_delete):
                new_size_b_buf += 1
                new_matrix_buf[index_a].append(new_matrix[index_a][index_b])

    new_matrix = new_matrix_buf
    new_size_b = int(new_size_b_buf/new_size_a)
    print("Iteration matrix = ",new_matrix)




A_size = 0
B_Size = 0
for row in new_matrix:
    A_size+=1
    buf = ""
    B_Size = 0
    for val in row:
        B_Size+=1
        buf+= str(val) + " "
    print(buf)

print("Result = = ",new_matrix)
print("New A size =",A_size)
print("New_B_Size =",B_Size)
P = []
Q = []
for i in range(0,A_size):
    P.append(Symbol('p' + str(i)))

for i in range(0,B_Size):
    Q.append(Symbol('q' + str(i)))

print(P)
print(Q)

Yy = Symbol('y')
P_system = [0]*A_size
Q_system = [0]*B_Size

for i in range(0,A_size):
    for j in range(0, B_Size):
        Q_system[i] += new_matrix[i][j]*P[j]
    Q_system[i] -= Yy
    #print("Q[",i,"]",Q_system[i])



for i in range(0,B_Size):
    for j in range(0, A_size):
        P_system[i] += new_matrix[j][i]*Q[j]
    P_system[i] -= Yy
    #print("P[",i,"]",P_system[i])

Q_system.append(-1)
P_system.append(-1)





Q.append(Yy)
P.append(Yy)
P_system,Q_system = Q_system,P_system
for j in range(0, B_Size):
    Q_system[len(Q_system) - 1] +=Q[j]
for j in range(0, A_size):
    P_system[len(P_system) - 1] +=P[j]

for i in range(0,A_size+1):
    print("Q[",i,"]",Q_system[i])
for i in range(0,B_Size+1):
    print("P[",i,"]",P_system[i])
print("Q_c =",Q)
print("P_c =",P)
print("Q_s =",Q_system)
print("P_s =",P_system)

print("\nResults:\n")

for solution_1 in linsolve(Q_system, Q):
    Result = solution_1
print("P = ",Result)

check_sym = 0
for j in range(0,A_size):
    print(P[j]," =",float(Result[j]))
    check_sym += float(Result[j])

print("Chek sym =",check_sym)
for solution_2 in linsolve(P_system, P):
    Result = solution_2
print("Q=",Result)

check_sym = 0
for j in range(0,A_size):
    print(Q[j]," =",float(Result[j]))
    check_sym+=float(Result[j])

print("Chek sym =",check_sym)
print("y =",float(Result[len(Result)-1]))

