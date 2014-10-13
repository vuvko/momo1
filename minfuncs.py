import numpy as np
import matplotlib.pyplot as plt


def plot_points(a, b, fa, fb):
    plt.plot(a, fa, 'rx')
    plt.plot(b, fb, 'go')
    plt.draw()
    #plt.show(block=0)


def min_parabolic_step(func, bounds, display=0):
    return bounds, num_oracle


def min_cubic_step(func, bounds, display=0):
    return bounds, num_oracle


def min_golden_step(func, bounds, display=0):
    k = (np.sqrt(5) - 1) / 2.0
    a, fa, b, fb, l, u, int_len = bounds
    new_int_len = k * int_len
    if display:
        print('  Old interval:', int_len)
        print('  New interval:', new_int_len)
        print('  Old points:', l, a, b, u)
    if fa > fb:
        l = a
        a = b
        b = l + new_int_len
        fa = fb
        fb = func(b)
    else:
        u = b
        b = a
        a = u - new_int_len
        fb = fa
        fa = func(a)
    if display:
        print('  New points:', l, a, b, u)
        plot_points(a, b, fa, fb)
    bounds = (a, fa, b, fb, l, u, new_int_len)
    return bounds, 1


def min_golden(func, bounds, tol_x=1e-5, max_iter=500, display=0):
    l, u = bounds
    k = (np.sqrt(5) - 1) / 2.0
    int_len = u - l
    new_int_len = k * int_len
    a = u - new_int_len
    fa = func(a)
    b = l + new_int_len
    fb = func(b)
    num_oracle = 2
    exitflag = 1
    for k in range(max_iter):
        if display:
            print('Iteration #', k)
        bounds = (a, fa, b, fb, l, u, int_len)
        bounds, num_oracle_plus = min_golden_step(func, bounds, display)
        a, fa, b, fb, l, u, int_len = bounds
        num_oracle = num_oracle + num_oracle_plus
        if int_len < tol_x:
            exit_flag = 0
            break
    if fa < fb:
        f_min = fa
        x_min = a
    else:
        fmin = fb
        x_min = b
    return x_min, f_min, exitflag, num_oracle


def wolfe_condition(f0, df0, a, dfa, c1, c2):
    flag = 0
    if fa <= f0 + c1 * a * df0:
        flag = flag + 1
    if np.abs(dfa) >= c2 * df0:
        flag = flag + 2
    return flag


def min_wolfe(func, alpha0, c1=1e-4, c2=0.1, alpha_max=100, 
  tol_interval=1e-8, max_iter=500, display=0):
    k = (np.sqrt(5) - 1) / 2.0
    a = 0
    f0, df0 = func(0)
    fa, dfa = f0, df0
    b = alpha0
    fb, dfb = func(b)
    x = w = v = (a + b) / 2.0
    fx, dfx = func(x)
    fw = fv = fx
    dfw = dfv = dfx
    exitflag = 2
    num_oracle = 3
    wolfe_flag0 = wolfe_condition(f0, df0, b, dfb, c1, c2)
    wolfe_flag1 = wolfe_condition(f0, df0, x, dfx, c1, c2)
    if wolfe_flag0 = 1 or wolfe_flag1 = 1:
        wolfe_flag = 1
    else:
        wolfe_flag = 0
    for k in range(max_iter):
        if display:
            print('Iteration #', k)
        if wolfe_flag: # ищем точку, где выполняется первое условие Вольфа
            
        else: # ищем точку, где вдобавок выполняется второе условие Вольфа
    return alpha_min, phi_min, exitflag, num_oracle
