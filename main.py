#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import check_grad

from minfuncs import min_golden, min_wolfe


def f1(x):
    b = 2
    y = - x / (x ** 2.0 + b)
    dy = (x ** 2 - b) / ((x ** 2 + b) ** 2.0)
    return y, dy


def f2(x):
    b = 0.01
    l = 39
    if x <= 1 - b:
        f0 = 1 - x
        df0 = -1
    elif x >= 1 + b:
        f0 = x - 1
        df0 = 1
    else:
        f0 = (0.5 / b) * (x - 1) ** 2 + b / 2
        df0 = (x - 1.0) / b
    y = f0 + ((2 - 2 * b) / (l * np.pi)) * np.sin(x * l * np.pi / 2)
    dy = df0 + (1 - b) * np.cos(x * l * np.pi / 2)
    return y, dy


def g(b):
    return np.sqrt(1 + b ** 2) - b


def f3(x, b1, b2):
    y = g(b1) * np.sqrt((1 - x) ** 2 + b2 ** 2) + g(b2) * np.sqrt(x ** 2 + b1 ** 2)
    dy = (g(b1) * (x - 1)) / np.sqrt((1 - x) ** 2 + b2 ** 2) + (g(b2) * x) / np.sqrt(x ** 2 + b1 ** 2)
    return y, dy


def f3_1(x):
    b1 = 0.001
    b2 = 0.001
    return f3(x, b1, b2)


def f3_2(x):
    b1 = 0.01
    b2 = 0.001
    return f3(x, b1, b2)


def f3_3(x):
    b1 = 0.001
    b2 = 0.01
    return f3(x, b1, b2)


def f4(x):
    y = -(x + 1) ** 2 + 1
    dy = -2 * x
    return y, dy


def f5(x):
    y = -np.sqrt(x + 1) + 1
    dy = -1.0 / (2 * np.sqrt(x + 1))
    return y, dy


def plot_func(func, a, b):
    x = np.arange(a, b, (b - a) / 100.0)
    Y = np.array([func(e) for e in x])
    print(Y.shape)
    plt.plot(x, Y[:, 0], 'b-')
    plt.draw()


def main():
    func = f1
    alpha0 = 0.5
    alpha_max = 10
    plt.figure()
    plot_func(func, 0, alpha_max)
    x_min, f_min, exitflag, num_oracle = min_wolfe(func, alpha0, c1=0.1, c2=0.21, alpha_max=alpha_max, display=1)
    print('Final result: x_min =', x_min, 'f_min =', f_min, 'exitflag =', exitflag, 'num_oracle =', num_oracle)
    plt.show()

if __name__ == '__main__':
    main()
