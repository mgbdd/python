import sympy as sp
import math
import argparse

E = 10**(-5)
D = 1

def define_left_border(inf, b, func):
    delta = b
    while func.subs(x, delta) >= 0:
        delta -= D
    return delta

def define_right_border(a, inf, func):
    delta = a
    while func.subs(x, delta) <= 0:

        delta += D
    return delta

def division_method(a, b, func):    
    if a == -math.inf:
        a = define_left_border(a, b, func)
    if b == math.inf:
        b = define_right_border(a, b, func)

    if func.subs(x, a) * func.subs(x, b) >= 0:
        raise ValueError("function has to change its sign on [a, b]")
    c = a
    while abs(b - a) > E:
        c = (a + b) / 2
        if func.subs(x, c) == 0:
            return c
        elif func.subs(x, a) * func.subs(x, c) < 0:
            b = c
        else:
            a = c
    return c

def negative_discriminant(roots, func): #функция монотонно возрастает, один корень
    r = 0
    if func.subs(x, 0) < -E: #f(0) < 0 -> r на (0, +oo)
        r = division_method(0, math.inf, func)
    elif func.subs(x, 0) > E: #f(0) > 0 -> r на (-oo, 0)
        r = division_method(-math.inf, 0, func)
    roots[float(r)] = 1

def zero_discriminant(roots, func, derivative): #три одинаковых корня
    x1 = sp.solve(derivative)[0]
    r = x1
    if func.subs(x, x1) > E:  #f(x1) > 0 -> r на (-oo, x1)
        r = division_method(-math.inf, x1, func)
    elif func.subs(x, x1) < -E: # f(x1) < 0 -> r на (x1, +oo)
        r = division_method(x1, math.inf, func)
    roots[r] = 3 

def positive_discriminant(roots, func, derivative): 
    x1, x2 = sp.solve(derivative)
    if x1 > x2:  #предполагаем, что x1 левее x2, x1 - max, x2 - min
        tmp = x1
        x1 = x2
        x2 = tmp
    
    #один корень кратности один
    if func.subs(x, x2) > E:   # f(x2) > 0 -> r на (-oo, x1)
        r = division_method(-math.inf, x1, func)
        roots[float(r)] = 1
    elif func.subs(x, x1) < -E: # f(x1) < 0 -> r на (x2, +oo)
        r = division_method(x2, math.inf, func)
        roots[float(r)] = 1

    #два различных корня, один из них кратности два
    elif abs(func.subs(x, x2)) < E and func.subs(x, x1) > E: # f(x2) = 0 & f(x1) > 0 -> r на (-oo, x1)
        roots[x2] = 2
        r = division_method(-math.inf, x1, func)
        roots[float(r)] = 1
    elif abs(func.subs(x, x1)) < E and func.subs(x, x2) < -E: # f(x1) = 0 & f(x2) < 0 -> r на (x2, +oo)
        roots[x1] = 2
        r = division_method(x2, math.inf, func)
        roots[float(r)] = 1

    #три различных корня
    elif func.subs(x, x1) > E and func.subs(x, x2) < -E: # f(x1) > 0 & f(x2) < 0 
        r1 = division_method(x1, x2, func)
        r2 = division_method(-math.inf, x1, func)
        r3 = division_method(x2, math.inf, func)
        roots[float(r1)] = 1
        roots[float(r2)] = 1
        roots[float(r3)] = 1 




if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Solve a cubic equation')
    parser.add_argument('a', type=float, help='The coefficient of x^3')
    parser.add_argument('b', type=float, help='The coefficient of x^2')
    parser.add_argument('c', type=float, help='The coefficient of x')
    parser.add_argument('d', type=float, help='Last coefficient')
    args = parser.parse_args()

    a, b, c, d = args.a, args.b, args.c, args.d

if a < 0:
    a*=(-1)
    b*=(-1)
    c*=(-1)
    d*=(-1)

x = sp.Symbol('x')
f = a*x**3 + b*x**2 + c*x + d
print("Your function: ", f)

derivative = sp.diff(f, x) 
print("Derivative of your function: ", derivative)

d = sp.discriminant(derivative, x)
print("Discriminant: ", f"{d:.3f}")

roots = dict()

if d < 0:
    negative_discriminant(roots, f)
elif d == 0.0:
    zero_discriminant(roots, f, derivative)
else:
    positive_discriminant(roots, f, derivative)

for key, value in roots.items():
    print(f"Root: {key:.3f}, multiplicity: {value}")






