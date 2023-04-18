from scipy.stats import bernoulli
from statistics import mean, variance

p = 0.5
n = 10
data_ans = []
data_bern = []
s = 100

for i in range(n):
    data_bern = bernoulli.rvs(size=s, p=p)
    print(data_bern)
    for j in range(s):
        if data_bern[j] == 1:
            data_ans.append(j+1)
            break

print(data_ans)

for el in data_ans:
    element_count = len([e for e in data_ans if e == el])
    print(el, element_count/len(data_ans))


print("Мат ожидание (теор):", 1/p)
print("Мат ожидание:", mean(data_ans))
print("Дисперсия (теор):", (1-p)/p**2)
print("Дисперсия:", variance(data_ans))