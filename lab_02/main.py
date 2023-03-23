import math
import datetime

from constants import Constants
from RungeKutta import RungeKutta
from tools import (
    Interpolation,
    Integral,
    draw_xy_graph,
    print_result_like_table
)

# Находим функции, используя интерполяцию на данные в тз таблицы
T0 = Interpolation(Constants.I, Constants.T0)
m = Interpolation(Constants.I, Constants.m)
sigma = Interpolation(Constants.T, Constants.Sigma)


def integrand_expression(x, I):
    """ Используется в функции Rp"""
    Tw = Constants.Tw
    return x * sigma(T0(I) + (Tw - T0(I)) * x ** m(I))


def Rp(I):
    lp = Constants.lp
    R = Constants.R
    PI = math.pi
    integral = Integral(integrand_expression, 0, 1)

    return lp / (2 * PI * R * R * integral(100, I))


def dIdt(t, U, I):
    return (U - (Constants.Rk + Rp(I)) * I) / Constants.Lk


def dUdt(t, U, I):
    return -I / Constants.Ck


def wrap_with_time(r_obj, order_accuracy):
    start = datetime.datetime.now()
    i_res, u_res = r_obj(order_accuracy)
    end = datetime.datetime.now()

    print(f'Порядок {order_accuracy} ---> OK')
    return i_res, u_res, end - start


def wrap_one_calculate(r_obj, order_accuracy):
    I, U, T = wrap_with_time(r_obj, order_accuracy)

    R = [(i[0], Rp(math.fabs(i[1]))) for i in I]
    IR = [(i[0], R[ind][1] * i[1]) for ind, i in enumerate(I)]

    # _T0 = [(i[0], T0(math.fabs(i[1]))) for i in I]
    # draw_xy_graph(_T0, color='orange', label=f'T0 ({order_accuracy})')

    draw_xy_graph(U, color='red', label=f'U ({order_accuracy})')
    draw_xy_graph(R, color='blue', label=f'R ({order_accuracy})')
    draw_xy_graph(IR, color='orange', label=f'IR ({order_accuracy})')
    draw_xy_graph(I, color='green', label=f'I ({order_accuracy})', show=True)

    return T


def main():
    r_obj = RungeKutta(dIdt, dUdt, 0, 0.5, 1400, 700e-6, 0.5e-6)

    t1 = wrap_one_calculate(r_obj, 1)
    t2 = wrap_one_calculate(r_obj, 2)
    t4 = wrap_one_calculate(r_obj, 4)

    print_result_like_table(
        '\nВРЕМЯ ВЫПОЛНЕНИЯ',
        [[t1 + t2 + t4, t1, t2, t4]],
        ['Всего', '1 порядок', '2 порядок', '4 порядок']
    )


if __name__ == '__main__':
    main()

