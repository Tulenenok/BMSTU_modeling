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


def integrand_expression(x, I):
    """ Используется в функции Rp"""

    # Находим функции, используя интерполяцию на данные в тз таблицы
    T0 = Interpolation(Constants.I, Constants.T0)
    m = Interpolation(Constants.I, Constants.m)
    sigma = Interpolation(Constants.T, Constants.Sigma)

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


def main():
    r_obj = RungeKutta(dIdt, dUdt, 0, 0.5, 1400, 700e-6, 0.5e-6)

    i1, u1, t1 = wrap_with_time(r_obj, 1)
    i2, u2, t2 = wrap_with_time(r_obj, 2)
    i4, u4, t4 = wrap_with_time(r_obj, 4)

    draw_xy_graph(u1, color='red', label='U (1)')
    draw_xy_graph(i1, color='green', label='I (1)', show=True)

    draw_xy_graph(u2, color='orange', label='U (2)')
    draw_xy_graph(i2, color='blue', label='I (2)', show=True)

    draw_xy_graph(u4, color='magenta', label='U (4)')
    draw_xy_graph(i4, color='purple', label='I (4)', show=True)

    print_result_like_table(
        '\nВРЕМЯ ВЫПОЛНЕНИЯ',
        [[t1 + t2 + t4, t1, t2, t4]],
        ['Всего', '1 порядок', '2 порядок', '4 порядок']
    )


if __name__ == '__main__':
    main()

