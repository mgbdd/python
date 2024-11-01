import spline_utils as su
from spline_utils import n
from spline_utils import l
from spline_utils import r

import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

def fabs(x):
    return np.abs(x)

def first_condition(spline_list, h, exp_list, const_list):

    x = sp.Symbol('x')

    for i in range(n):
        const_list.append(fabs(l + h * i))
        const_list.append(fabs(l + h * (i + 1)))

        s = spline_list[i].subs(x, l + h * i)
        s = sp.sympify(s).evalf(3)
        exp_list.append(s)

        s = spline_list[i].subs(x, l + h * (i+1))
        s = sp.sympify(s).evalf(3)
        exp_list.append(s)

    coef_matrix = su.create_null_matrix()

    for i in range(2*n):
        su.get_coefs(exp_list[i], i, coef_matrix)
    
    return coef_matrix

def second_condition(spline_list, h, coef_matr, const_list):
    x = sp.Symbol('x')
    first_deriatives = list()
    second_deriatives = list()
    for spline in spline_list:
        first_deriatives.append(sp.diff(spline, x))
        second_deriatives.append(sp.diff(first_deriatives[-1], x))

    for i in range(n - 1):
        ex1 = first_deriatives[i].subs(x, l + (i+1)*h) - first_deriatives[i+1].subs(x, l + (i+1)*h)
        ex1 = sp.sympify(ex1).evalf(3)
        ex2 = second_deriatives[i].subs(x, l + (i+1)*h) - second_deriatives[i+1].subs(x, l + (i+1)*h)
        ex2 = sp.sympify(ex2).evalf(3)
        su.get_coefs(ex1, 2*n + 2*i, coef_matr)
        su.get_coefs(ex2, 2*n + 2*i + 1, coef_matr)
    
    for i in range(2*n, 4*n - 2):
        const_list.append(0)

def third_condition(spline_list, coef_matr, const_list):
    const_list.append(0)
    const_list.append(0)

    x = sp.Symbol('x')
    s0_1 = sp.diff(spline_list[0], x)
    s0_2 = sp.diff(s0_1, x)
    su.get_coefs(s0_2.subs(x, l), 4*n - 2, coef_matr)
    #su.get_coefs(spline_list[0].subs(x, l), 4*n - 2, coef_matr)
    #ex1 = spline_list[0].subs(x, l) - spline_list[n-1].subs(x, r)
    #su.get_coefs(ex1, 4*n - 2, coef_matr)

    s0_1 = sp.diff(spline_list[n-1], x)
    s0_2 = sp.diff(s0_1, x)
    su.get_coefs(s0_2.subs(x, r), 4*n - 1, coef_matr)
    #su.get_coefs(spline_list[n-1].subs(x, r), 4*n - 1, coef_matr)
    #ex1 = sp.diff(spline_list[0].subs(x, l), x) - sp.diff(spline_list[n-1].subs(x, r), x)
    #su.get_coefs(ex1, 4*n - 2, coef_matr)

def solve_system(coef_matr, const_list, h):
    solution = np.linalg.solve(coef_matr, const_list)
    #print(solution)

    plt.figure(figsize=(10, 6)) 
    for i in range(n):
        a = solution[i * 4]
        b = solution[i * 4 + 1]
        c = solution[i * 4 + 2]
        d = solution[i * 4 + 3]
        x = np.linspace(l + i * h, l + (i + 1) * h)
        curr = l + i * h
        y = a * (x - curr) ** 3 + b * (x - curr) ** 2 + c * (x - curr) + d
        plt.plot(x, y, label=f'Spline {i + 1}')

    plt.title('Кусочная функция из кубических сплайнов')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axhline(0, color='black', lw=0.5)
    plt.axvline(0, color='black', lw=0.5)
    plt.grid()
    plt.savefig('spline.png')



def main():

    x = np.linspace(l, r, 100)
    y = fabs(x)
    su.draw_graph(x, y, "orig", "f(x) = |x|")

    h = (abs(l) + abs(r)) / n



    spline_list = su.calculate_splines(h)
    exp_list = list() # список с уравнениями, в которых неизвестные - это коэффициенты сплайнов
    const_list = list() # список со свободными членами будущей матрицы коэффициентов

    print("SPLINE LIST")
    for i in range(n):
        print(spline_list[i])
 
    coef_matr = first_condition(spline_list, h, exp_list, const_list) #получаем матрицу с заполненными 2n строками 

    print("EXP LIST")
    for i in range(2*n):
        print(exp_list[i])

    #for i in range(2*n):
        #print(coef_matr[i])
#second_condition(spline_list, h, coef_matr, const_list) # 2(n-1) строчек
#third_condition(spline_list, coef_matr, const_list) # последние 2 строчки 

    #for i in range(4*n):
        #print(coef_matr[i])

    solve_system(coef_matr, const_list, h)


#в третьем условии взяла =0, но могут быть и другие случаи, просто этот самый простой



if __name__ == "__main__":
    main()


# сколько интервалов -> столько сплайнов (n)
# формула сплайна: S_i = a_i + b_i * (x - x_j) + c_i * (x - x_j)^2 + d_i * (x - x_j)^3, i - номер сплайна  
# условия сплайнов:
# 1. проходит через узловые точки -> подставить точки -> получить n уравнений с 4n неизвестными коэффициентами
# 2. в стыках сплайнов должна обеспечиваться гладкость -> для каждого сплайна находим 1 и 2 производные -> подставляем точки стыка (S'_0 = S'_1, S"_0 = S"_1)
# 3.проверка краевых точек
# получим 4n уравнений -> найдем все неизвестные
# 4. построить графики

