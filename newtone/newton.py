import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

l = -1
r = 1
n = 10 # количество интервалов
# 10 интервалов -> 11 точек
# точки: x0, x1, ... x10

def calc_delta(h, k):
    sum = 0
    for i in range(k + 1):
        mult = 1
        for j in range(k + 1):
            if i == j:
                continue
            else: 
                mult *= (l + i * h) - (l + j * h) #(x_i - x_j)
        sum += abs(l + i * h) / mult
    
    return sum
    
def build_polynomial(h):
    x = sp.Symbol("x")
    f = abs(l)
    for k in range(1, n + 1):
        coef = calc_delta(h, k)
        big_mult = 1
        for i in range(k):
            mult = (x - (l + i * h))
            big_mult *= mult
        f += coef * big_mult

    f = sp.expand(f)
    print(f)
    return sp.lambdify(x, f, 'numpy')

def draw_func(polynomial):
    x = np.linspace(1*l, 1*r, n+1)
    y = polynomial(x)
    plt.figure(figsize=(16, 10))  
    plt.plot(x, y, color='blue')  
    plt.title('График многочлена в форме Ньютона')  
    plt.xlabel('x')  
    plt.ylabel('y')  
    plt.axhline(0, color='black',linewidth=0.5, ls='--')  # Горизонтальная линия y=0
    plt.axvline(0, color='black',linewidth=0.5, ls='--')  # Вертикальная линия x=0
    plt.grid(color = 'gray', linestyle = '--', linewidth = 0.5)  # Сетка
    plt.savefig('newton.png')
    

def main():
    h = (abs(l) + abs(r)) / n

    x = sp.Symbol("x")
    f = str((abs(l)))
    coef = -1
    mult = str((x - l))
    f += str(coef) + "*" + "(" + mult + ")"
    #print(f)
    f = sp.sympify(f)
    #print(f.subs(x, 5))


    polynomial = build_polynomial(h)
    draw_func(polynomial)

if __name__ == "__main__":
    main()