print("Data max/min = ")
Aray_data = []
Aray_data = input().split()

print("Functions = ", end=' ')
Count_Functions = int(input())

array = []
Buf_array = []
for i in range(Count_Functions):
    array.append([])
    print("Input krits for ", i+1, "rhtow space bar = ")
    Buf_array = input().split()
    for Element in Buf_array:
        array[i].append(float(Element))


Max_Min_vals_array = []
for Krit in Aray_data:
    if Krit == "min":
        Max_Min_vals_array.append(1000)
    else:
        Max_Min_vals_array.append(-1000)

for i in range(Count_Functions):
    j = 0
    for Krit in Aray_data:
        if Krit == "min":
            if Max_Min_vals_array[j] > array[i][j]:
                Max_Min_vals_array[j] = array[i][j]
        if Krit == "max":
            if Max_Min_vals_array[j] < array[i][j]:
                Max_Min_vals_array[j] = array[i][j]
        j + = 1

Result = []
Check = False
for i in range(Count_Functions):
    Check = False
    j = 0
    for Krit in Aray_data:
        if array[i][j] == Max_Min_vals_array[j]:
            Check = True
    if Check:
        Result.append(array[i])
    j + = 1

print("Paretto = ")
for res in Result:
    print(res)
