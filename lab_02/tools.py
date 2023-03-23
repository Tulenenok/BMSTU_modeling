import matplotlib.pyplot as plt
import math
from scipy.interpolate import interp1d
from scipy.integrate import quad
from prettytable import PrettyTable
from constants import Constants


class Interpolation:
    def __init__(self, x_list, y_list):
        self.x = x_list
        self.y = y_list

        # self.f = interp1d(self.x, self.y)

    # def __call__(self, x):
    #     return self.f(x)

    def __call__(self, x0):
        x_log = [math.log(x) for x in self.x]
        y_log = [math.log(y) for y in self.y]

        x0 = math.log(x0)
        y0 = 0

        for i in range(len(x_log)):
            l = 1
            for j in range(len(x_log)):
                if i != j:
                    l *= (x0 - x_log[j]) / (x_log[i] - x_log[j])
            y0 += y_log[i] * l

        return math.exp(y0)


class Integral:
    def __init__(self, f, a, b):
        self.a = a
        self.b = b
        self.f = f

    # def var2(self, n, I):
    #     return quad(self.f, self.a, self.b, kind='cubic')

    def __call__(self, n, I):
        w = (self.b - self.a) / n
        integral = 0

        for i in range(n):
            x1 = self.a + w * i
            x2 = self.a + w * (i + 1)

            integral += 0.5 * (x2 - x1) * (self.f(x1, I) + self.f(x2, I))

        return integral


def generate_grid_use_step(start, end, step):
    grid = []
    while start <= end:
        grid.append(start)
        start += step

    return grid


def generate_grid_use_cnt(start, end, n):
    step = (end - start) / n
    return generate_grid_use_step(start, end, step)


def draw_graph(x, y, color='b', label='', show_points=False, show=False):
    plt.plot(x, y, color, label=label)

    if show_points:
        plt.plot(x, y, color + 'o')

    ax = plt.gca()
    # ax.spines['left'].set_position('center')
    # ax.spines['bottom'].set_position('center')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    plt.grid()
    plt.legend()

    if show:
        plt.show()


def draw_xy_graph(xy, color='b', label='', show_points=False, show=False):
    x = [i[0] for i in xy]
    y = [i[1] for i in xy]

    draw_graph(x, y, color, label, show_points, show)


def print_result_like_table(label, lst, schema):
    print(label)

    tb = PrettyTable()
    tb.field_names = schema

    for dt in lst:
        tb.add_row(dt)

    print(tb, '\n')


def test_interpolation():
    """ Если понадобиться писать свою интерполяцию, можно затестить ее этом кодом,
        сверив с графиками их папки img
    """
    T0 = Interpolation(Constants.I, Constants.T0)
    m = Interpolation(Constants.I, Constants.m)
    sigma = Interpolation(Constants.T, Constants.Sigma)

    I = generate_grid_use_step(0.5, 1200, 0.01)
    T = generate_grid_use_step(4000, 14000, 0.001)

    T0_list = [T0(i) for i in I]
    m_list = [m(i) for i in I]

    draw_graph(I, T0_list, label='T0', show=True)
    draw_graph(I, m_list, color='m', label='m', show=True)
    draw_graph(T, sigma(T), color='y', label='sigma', show=True)
