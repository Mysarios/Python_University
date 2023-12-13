import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import rc

Stranght = 1.4
high = 100.0

coef_turb_z = 5.0
k0 = 0.5
u = 5.0
coef_turb_y = k0 * u

x, y = np.mgrid[0.1:10:100j, -5:5:100j]
z = 100

result = (Stranght / (4 * math.pi * x * math.sqrt(coef_turb_y * coef_turb_z))
          * np.exp(-u * y ** 2 / (4 * coef_turb_y * x))
          * (np.exp(-u * (z + high) ** 2 / (4 * coef_turb_z * x)) + np.exp(-u * (z - high) ** 2 / (4 * coef_turb_z * x)))
          )

# расчет для y = 0; вертикальная плоскость
x1, z1 = np.mgrid[0.1:10:100j, 90:110:100j]
y1 = 0.
result_y0 = (Stranght / (4 * math.pi * x1 * math.sqrt(coef_turb_y * coef_turb_z))
             * np.exp(-u * y1 ** 2 / (4 * coef_turb_y * x1))
             * (np.exp(-u * (z1 + high) ** 2 / (4 * coef_turb_z * x1)) + np.exp(-u * (z1 - high) ** 2 / (4 * coef_turb_z * x1)))
             )

# определение линий уровня и цветов
Line_highs= [0.002, 0.004, 0.006, 0.008, 0.010, 0.012, 0.016, 0.020, 0.030, 0.040, 0.300]
Color_lines= ['cyan', 'turquoise', 'teal', 'steelblue', 'cornflowerblue', 'mediumslateblue',
       'indigo', 'darkmagenta', 'orchid', 'lightpink', 'hotpink']

# выводим два графика рядом

plt.contourf(x, y, result, levels = Line_highs, colors=Color_lines)
plt.xlabel('$x$, м')
plt.ylabel('$y$, м')
plt.title('$а$) $z$ = 100 м', loc='left')
plt.show()


plt.contourf(x1, z1, result_y0, levels = Line_highs, colors=Color_lines)
plt.xlabel('$x$, м')
plt.ylabel('$z$, м')
plt.title('$б$) $y$ = 0', loc='left')

plt.show()
#plt.savefig('d:\Lea\Ecolog\dif_analit.png')