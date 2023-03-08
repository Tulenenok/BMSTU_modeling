import matplotlib.pyplot as plt
from math import e
from scipy.integrate import quad
from prettytable import PrettyTable


def derivative(x, u):
    """ Производная функции (в данном случае dx/du = u^2 + x) """
    return u ** 3 + 2 * x * u


def analytical_function(y):
    """ Аналитическое решение уравнения """
    return e ** (y ** 2) - (y ** 2) / 2 - 0.5


def find_solution_for_equation(func, x0):
    """ Найти решение уравнения 0 = func"""
    x = x0
    h = 0.00001

    df = (func(x + h) - func(x)) / h
    for i in range(1000):
        x = x - func(x) / df
    return x


def solve_euler_2(der, x0, y0, hy, y_end):
    result = []
    x, y = x0, y0

    while hy > 0 and y <= y_end or hy < 0 and y >= y_end:
        result.append((x, y))

        tmp = lambda xi: xi - hy * der(xi, y + hy) - x
        x = find_solution_for_equation(tmp, x)
        y += hy

    return result


def solve_euler(der, x0, y0, hy, y_end):
    """
    Явный метод Эйлера

    :param der: производная функции dx/du или du/dx
    :param x0: координата x, с которой хотим начать расчет
    :param y0: координата y, с которой хотим начать расчет
    :param hx: шаг по оси x
    :param x_end: координата x, на которой хотим закончить

    :return: массив точек вида (x, y)
    """

    result = []
    x, y = x0, y0

    while hy > 0 and y <= y_end or hy < 0 and y >= y_end:
        result.append((x, y))
        x += hy * der(x, y)
        y += hy

    return result


def solve_analytical(func, x0, y0, hy, y_end):
    """
    Аналитическое решение (получаем точки для построения графика)

    :param func: решение уравнения.
    :param x0: координата x, с которой хотим начать расчет.
    :param y0: координата y, с которой хотим начать расчет.
    :param hy: шаг по оси y.
    :param y_end: координата y, на которой хотим закончить.

    :return: массив точек вида (x, y)
    """

    result = []
    x = x0
    y = y0

    while (hy > 0 and y <= y_end) or (hy < 0 and y >= y_end):
        result.append((x, y))
        y += hy
        x = func(y)

    return result


def solve_picar(der, x0, y0, hy, y_end, approx):

    def get_integral(_der, xi, si, nu=0, ne=0):
        if si == 0:
            return nu

        tmp = lambda param: _der(param, get_integral(_der, param, si - 1, nu, ne))
        integral = quad(tmp, ne, xi)

        return nu + integral[0]

    result = []
    x, y = x0, y0

    while hy > 0 and y <= y_end or hy < 0 and y >= y_end:
        result.append((x, y))
        x = get_integral(der, y, approx, nu=x0, ne=y0)
        y += hy

    return result


def draw_graph(xy, color='b', label='', show_points=False, need_sort=True, show=False):
    """
    Рисует граф по списку точек

    :param xy: массив точек вида (x, y)
    :param color: цвет линии
    :param label: легенда
    :param show_points: нужно ли отдельно отмечать точки (по умолчанию -- нет)
    :param need_sort: нужна ли сортировка точек по оси x (по умолчанию -- да)
    :param show: вывести сразу (по умолчанию нет)

    :return :
    """

    if need_sort:
        xy.sort(key=lambda i: i[0])

    x = [i[0] for i in xy]
    y = [i[1] for i in xy]

    plt.plot(x, y, color, label=label)

    if show_points:
        plt.plot(x, y, color + 'o')

    ax = plt.gca()
    # ax.spines['left'].set_position('center')
    ax.spines['bottom'].set_position('center')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()

    if show:
        plt.show()


plt.xlim([0, 10])
plt.ylim([-3, 3])

X = 0.5
Y = 0
STEP = 0.1
END = 3

t = PrettyTable()

""" Аналитическое """
analytical_points = solve_analytical(analytical_function, X, Y, -STEP, -END)
draw_graph(analytical_points, color='m', label='analytical')

analytical_points = solve_analytical(analytical_function, X, Y, STEP, END)
draw_graph(analytical_points, color='m')

""" Эйлер """
euler_points = solve_euler(derivative, X, Y, -STEP, -END)
draw_graph(euler_points, label='explicit Euler')

euler_points = solve_euler(derivative, X, Y, STEP, END)
draw_graph(euler_points)

""" Неявный Эйлер"""
implicit_points = solve_euler_2(derivative, X, Y, STEP, END)
draw_graph(implicit_points, color='g', label='implicit Euler')


""" Заполняем таблицу """
t.add_column("X", [i[0] for i in euler_points])
t.add_column("Analytical", [i[1] for i in analytical_points[:len(euler_points)]])
t.add_column("Euler", [i[1] for i in euler_points])


""" Пикар """
colors = ['#F9ED69', '#F08A5D', '#B83B5E', '#6A2C70']
for i in range(1, 5):
    # picar_points = solve_picar(derivative, X, Y, -STEP, -END, i)
    picar_points = solve_picar(derivative, X, Y, STEP, END, i)
    t.add_column(f"Piсard {i}", [i[1] for i in picar_points])
    draw_graph(picar_points, color=colors[i - 1], label=f'Picard {i}')

print(t)
plt.show()

