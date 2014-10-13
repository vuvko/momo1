import numpy as np
import matplotlib.pyplot as plt

from minfuncs import min_golden, min_wolfe


def f1(x):
    return 0.2 * x * np.log(x) + (x - 2.3) ** 2


def plot_func(func, a, b):
    x = np.arange(a, b, (b - a) / 100.0)
    y = func(x)
    plt.plot(x, y, 'b-')
    plt.draw()


def main():
    l = 0.5
    u = 2.5
    plt.figure()
    plot_func(f1, l, u)
    x_min, f_min, exitflag, num_oracle = min_golden(f1, (l, u), display=1)
    print('Final result: x_min =', x_min, 'f_min =', f_min, 'num_oracle =', num_oracle)
    plt.show()

if __name__ == '__main__':
    main()
