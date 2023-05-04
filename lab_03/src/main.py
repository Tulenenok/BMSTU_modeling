import matplotlib.pyplot as plt

from solver import (
    solve,
    Input,
)
from tools import TableFunc


def draw_graph(x, y, color='b', label='', show_points=False, show=False):
    plt.plot(x, y, color, label=label)

    if show_points:
        plt.plot(x, y, color + 'o')

    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.grid()
    plt.legend()

    if show:
        plt.show()


def main():
    lambda_t_list = [300, 500, 800, 1100, 2000, 2400]
    lambda_list = [1.36 * 10**-2, 1.63 * 10**-2, 1.81 * 10**-2, 1.98 * 10**-2, 2.5 * 10**-2, 2.74 * 10**-2]

    k_t_list = [293, 1278, 1528, 1677, 2000, 2400]
    k_list = [2.0 * 10**-2, 5.0 * 10**-2, 7.8 * 10**-2, 1.0 * 10**-1, 1.3 * 10**-1, 2.0 * 10**-1]

    handle = Input(TableFunc(lambda_t_list, lambda_list), TableFunc(k_t_list, k_list))

    handle.n = 1.4
    handle.r0 = 0.35
    handle.r_max = 0.5
    handle.steps = 10000
    handle.t0 = 300
    handle.sigma = 5.668e-12
    handle.f_0 = 100
    handle.alpha = 0.05
    handle.eps1 = 10 ** -4
    handle.eps2 = 10 ** -4

    result = solve(handle)

    print('F1 = ', result.f1)
    print('F2 = ', result.f2)
    print('Iter = ', result.iter)

    draw_graph(result.radius, result.temperature, 'm', label='Зависимость Т(r)', show=True)


if __name__ == "__main__":
    main()
