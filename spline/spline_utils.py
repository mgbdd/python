import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

n = 10
l = -1
r = 1

def calculate_equation(x, h, ind):
    i = str(ind)
    x0 = l + h * ind
    f = sp.Symbol('a' + i)*(x - x0)**3 + sp.Symbol('b' + i)*(x - x0)**2 + sp.Symbol('c' + i)*(x - x0) + sp.Symbol('d' + i)
    f = sp.expand(f)
    rounded_expr = sp.sympify(f).evalf(3)
    return rounded_expr

def calculate_splines(h):
    x = sp.Symbol('x')
    spline_list = list()
    for i in range(n):
        spline_list.append(calculate_equation(x, h, i))

    return spline_list

def draw_graph(x, y, name, func):
    plt.figure(figsize=(8, 5))  
    plt.plot(x, y, color='blue')  
    plt.title('График функции ' + func)  
    plt.xlabel('x')  
    plt.ylabel('y')  
    plt.axhline(0, color='black',linewidth=0.5, ls='--')  # Горизонтальная линия y=0
    plt.axvline(0, color='black',linewidth=0.5, ls='--')  # Вертикальная линия x=0
    plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)  # Сетка
    plt.savefig(name + '.png')

def create_null_matrix():
    #matrix = [[0 for _ in range(4*n)] for _ in range(4*n)]
    matrix = np.zeros((4*n, 4*n))
    return matrix
              
def insert_coef(coef_matrix, exp_ind, spline_ind, letter, coef): # порядковый номер сплайна
    pos = None
    if letter == 'a':
        pos = 0
    elif letter == 'b':
        pos = 1
    elif letter == 'c':
        pos = 2
    else:
        pos = 3  
    pos += spline_ind * 4
    coef_matrix[exp_ind][pos] = coef

def get_coefs(exp, exp_ind, coef_matrix): # ind - порядковый номер сплайна 
    f = str(exp)
    start = 0
    next_coef_neg_flag = False
    for i in range(len(f)):
        if ord(f[start]) >= 97 and ord(f[start]) <= 100:
            insert_coef(coef_matrix, exp_ind, int(f[start + 1]), f[start], 1)
            if start + 5 <= len(f):
                start += 5
            else:
                break

        if f[i] == '*':
            if next_coef_neg_flag == False:
                insert_coef(coef_matrix, exp_ind, int(f[i + 2]), f[i + 1], float(f[start:i]))
            else:
                insert_coef(coef_matrix, exp_ind, int(f[i + 2]), f[i + 1], (-1)*float(f[start:i]))
                next_coef_neg_flag = False
        if f[i] == 'e':
            s_ind = None
            for j in range(i, len(f)):
                if f[j] == '*':
                    s_ind = j + 2
                    break
            insert_coef(coef_matrix, exp_ind, int(f[s_ind]), f[i + 5], float(f[start:(s_ind-2) ]))
        if f[i] == '+' :
            start = i + 2
        if f[i] == '-' and i != 0 and f[i - 1] != 'e':
            next_coef_neg_flag = True
            start = i + 2

     





# посмотреть после получения результата, стоит ли коэффициент с числом е заменять на ноль


#при n>10 вырожденная матрица коэффициентов