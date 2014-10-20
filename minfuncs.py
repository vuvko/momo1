#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt


def plot_points(a, b, fa, fb):
    plt.plot(a, fa, 'rx')
    #plt.plot(b, fb, 'go')
    plt.draw()
    #plt.show(block=0)


def plot_wolfe(f0, df0, alpha_max, c1):
    x = np.array([0, alpha_max])
    y = np.array([f0, f0 + c1 * alpha_max * df0])
    plt.plot(x, y, 'r-')
    plt.draw()


def min_parabolic_step_3p(func, bounds, display=0):
    x1, fx1, x2, fx2, x3, fx3 = bounds
    A = np.zeros((3, 3))
    A[0, :] = np.array([x1 ** 2, x1, 1])
    A[1, :] = np.array([x2 ** 2, x2, 1])
    A[2, :] = np.array([x3 ** 2, x3, 1])
    b = np.array([fx1, fx2, fx3])
    try:
        coeff = np.linalg.solve(A, b)
    except np.linalg.linalg.LinAlgError:
        coeff = np.ones(3)
        coeff[1] = -(x1 + x2)
    x_min = -coeff[1] / (2.0 * coeff[0])
    f_min = func(x_min)
    return x_min, f_min, 1


def min_parabolic_step_2p(func, bounds, display=0):
    x1, fx1, dfx1, x2, fx2, dfx2 = bounds
    A = np.zeros((3, 3))
    A[0, :] = np.array([x1 ** 2, x1, 1])
    A[1, :] = np.array([x2 ** 2, x2, 1])
    A[2, :] = np.array([2 * x1, 1, 0])
    b = np.array([fx1, fx2, dfx1])
    try:
        coeff = np.linalg.solve(A, b)
    except np.linalg.linalg.LinAlgError:
        coeff = np.ones(3)
        coeff[1] = -(x1 + x2)
    x_min = -coeff[1] / (2.0 * coeff[0])
    f_min = func(x_min)
    return x_min, f_min, 1


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


def wolfe_condition(f0, df0, a, fa, dfa, c1, c2):
    flag = 0
    if fa <= f0 + c1 * a * df0:
        flag = flag + 1
    if np.abs(dfa) <= c2 * abs(df0):
        flag = flag + 2
    return flag


def min_wolfe(func, alpha0, c1=1e-4, c2=0.1, alpha_max=100, 
  tol_interval=1e-8, max_iter=500, display=0):
    eps = tol_interval
    k = (np.sqrt(5) - 1) / 2.0
    a = 0
    f0, df0 = func(0)
    fa, dfa = f0, df0
    b = alpha0
    fb, dfb = func(b)
    x = v = w = b / 2.0
    fw, dfw = func(w)
    fx = fv = fw
    dfx = dfv = dfw
    exitflag = 2
    num_oracle = 3
    wolfe_flag = wolfe_condition(f0, df0, b, fb, dfb, c1, c2)
    if wolfe_flag == 1 and dfb < 0:
        wolfe_flag = -1
    alpha_min = b
    phi_min = fb
    
    plot_wolfe(f0, df0, alpha_max, c1)
    
    int_len = b - a
    for k in range(max_iter):
        if display:
            print('Iteration #', k)
            print('  a =', a, 'b =', b, 'w =', w, 'v =', v)
            print('  x =', x, 'fx =', fx, 'dfx =', dfx)
        if wolfe_flag == -1:
            if display:
                print('  Going right')
            # ищём правый край
            x = b * 2
            if x > alpha_max:
                alpha_min = b
                phi_min = fb
                exitflag = 3
                break
            fx, dfx = func(x)
            num_oracle = num_oracle + 1
            wolfe_flag = wolfe_condition(f0, df0, x, fx, dfx, c1, c2)
            int_len = x - b
            if wolfe_flag == 1 and dfx < 0 and fx < fb:
                wolfe_flag = -1
                b, fb, dfb = x, fx, dfx
            else:
                t, ft, dft = b, fb, dfb
                b, fb, dfb = x, fx, dfx
                x, fx, dfx = t, ft, dft
        else:
            # ищем точку, где выполняются условия Вольфа
            u1_flag = 0
            u2_flag = 0
            if x != v and dfx != dfv:
                if display:
                    print('  u1:')
                # параболическая аппроксимация
                u1, fu1, no_plus = min_parabolic_step_2p(func, (x, fx, dfx, v, fv, dfv), display)
                num_oracle = num_oracle + 1
                if (u1 > a + eps and u1 < b - eps and
                  abs(u1 - x) < (int_len / 1.5) and
                  not (wolf_flag != 1 and u1 > x)):
                    if display:
                        print('     yep.')
                    u1_flag = 1
            if x != w and dfx != dfw:
                if display:
                    print('  u2:')
                # параболическая аппроксимация
                u2, fu2, no_plus = min_parabolic_step_2p(func, (x, fx, dfx, w, fw, dfw), display)
                num_oracle = num_oracle + 1
                if (u2 > a + eps and u2 < b - eps and
                  abs(u2 - x) < (int_len / 1.5) and
                  not (wolfe_flag != 1 and u2 > x)):
                    if display:
                        print('     yep.')
                    u2_flag = 1
            if u1_flag and u2_flag:
                if abs(u2 - x) < abs(u1 - x):
                    if display:
                        print('  get u2')
                    u = u2
                else:
                    if display:
                        print('  get u1')
                    u = u1
            elif u1_flag:
                if display:
                    print('  get u1')
                u = u1
            elif u2_flag:
                if display:
                    print('  get u2')
                u = u2
            else:
                if display:
                    print('  No u:')
                if dfx > 0 or wolfe_flag != 1:
                    if display:
                        print('       [a, x]')
                    u = (a + x) / 2.0
                else:
                    if display:
                        print('       [x, b]')
                    u = (x + b) / 2.0
            if abs(u - x) < eps:
                u = x + np.sign(u - x) * eps
            int_len = abs(u - x)
            fu, dfu = func(u)
            if display:
                print('  u =', u, 'f_u =', fu)
            num_oracle = num_oracle + 1
            wolfe_flag_u = wolfe_condition(f0, df0, u, fu, dfu, c1, c2)
            if fu < fx:
                if u < x:
                    b = x
                elif wolfe_flag == 1:
                    a = x
                x, fx, dfx = u, fu, dfu
            else:
                if u < x:
                    if wolfe_flag_u == 1:
                        a = u
                    else:
                        x, fx, dfx = u, fu, dfu
                else:
                    b = u
            wolfe_flag = wolfe_condition(f0, df0, x, fx, dfx, c1, c2)
        if display:
            print('  wolfe:', wolfe_flag)
        if wolfe_flag == 3:
            alpha_min = x
            phi_min = fx
            exitflag = 0
            break
        if int_len < 2 * eps:
            alpha_min = x
            phi_min = fx
            exitflag = 1
            break
        if v != x:
            w, fw, dfw = v, fv, dfv
            v, fv, dfv = x, fx, dfx
        if display:
            plot_points(x, 0, fx, 0)
    if display:
        print('Done.')
        plt.plot(alpha_min, phi_min, 'go')
        plt.draw()
    return alpha_min, phi_min, exitflag, num_oracle
