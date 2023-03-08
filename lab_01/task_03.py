import matplotlib.pyplot as plt
from math import e
from scipy.integrate import quad
from prettytable import PrettyTable


def derivative(x, u):
    """ Производная функции """
    return x ** 2 + u ** 2


def solve_euler(der, x0, y0, hx, x_end):
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

    while hx > 0 and x <= x_end or hx < 0 and x >= x_end:
        result.append((x, y))
        y += hx * der(x, y)
        x += hx

    return result


def solve_picar(der, x0, y0, hx, x_end, approx):

    def get_integral(_der, xi, si, nu=0, ne=0):
        if si == 0:
            return nu

        tmp = lambda param: _der(param, get_integral(_der, param, si - 1, nu, ne))
        integral = quad(tmp, ne, xi)

        return nu + integral[0]

    result = []
    x, y = x0, y0

    while hx > 0 and x <= x_end or hx < 0 and x >= x_end:
        result.append((x, y))
        y = get_integral(der, x, approx, nu=y0, ne=x0)
        x += hx

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
    # ax.spines['bottom'].set_position('center')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.legend()

    if show:
        plt.show()


plt.xlim([0, 1.5])
plt.ylim([0, 1.5])

X = 0
Y = 0
STEP = 0.1
END = 3

t = PrettyTable()

""" Эйлер """
# euler_points = solve_euler(derivative, X, Y, -STEP, -END)
euler_points = solve_euler(derivative, X, Y, STEP, END)
draw_graph(euler_points, label='explicit Euler')

t.add_column("X", [i[0] for i in euler_points])
t.add_column("Euler", [i[1] for i in euler_points])

colors = ['#F9ED69', '#F08A5D', '#B83B5E', '#6A2C70']
for i in range(1, 5):
    # picar_points = solve_picar(derivative, 1, 0, -0.1, -3, i)
    picar_points = solve_picar(derivative, X, Y, STEP, END, i)
    t.add_column(f"Piсard {i}", [i[1] for i in picar_points])
    draw_graph(picar_points, color=colors[i - 1], label=f'Picard {i}')

print(t)
plt.show()



