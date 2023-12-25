import multiprocessing
from threading import Thread
from symengine import *
from sympy import *
import sympy
from scipy.integrate import trapz
import numpy as np
from sympy import factor_terms
from scipy import integrate
import time

Yy = Symbol('y')
Xx = Symbol('x')
Tt = symbols('t')
low = 0
high = 1
n = 5000
fun = []
quadratic_expression_1 = 5*Xx + 5*Yy
quadratic_expression_2 = 22*Xx**10 + 5*Yy
fun.append(quadratic_expression_1)
fun.append(quadratic_expression_2)
def func(f):
    print("f =",f)
    queue = multiprocessing.Queue()
    timer = time.time()
    p1 = multiprocessing.Process(target=intagrete, args=(f, 0,queue))
    p2 = multiprocessing.Process(target=intagrete, args=(f, 1,queue))
    p1.start()
    p2.start()
    p1.join()
    p2.join()
    print('Time %.6f' % (time.time() - timer))
    while queue.qsize() > 0:
        print(queue.get())
def intagrete_Es(fun,queue,type):
    Result = sympy.integrate(fun, (Xx, 0, 10))
    fun[i] = Result
    queue.put(fun)
    print("result =",Result)

if __name__ == '__main__':
    a = 1