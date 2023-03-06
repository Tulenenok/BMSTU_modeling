import matplotlib.pyplot as plt
from math import e

plt.xlim([-5, 5])
plt.ylim([-5, 5])


class Graph:
    @staticmethod
    def show_list(x_y, color='b', label='', show_points=False):
        """ Вывести граф по списку точек вида (x, y)"""
        x_y.sort(key=lambda i: i[0])

        x = [i[0] for i in x_y]
        y = [i[1] for i in x_y]

        plt.plot(x, y, color, label=label)

        if show_points:
            plt.plot(x, y, color + 'o')

        ax = plt.gca()
        ax.spines['left'].set_position('center')
        ax.spines['bottom'].set_position('center')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.legend()

    @staticmethod
    def draw_all():
        plt.show()


class ODY:
    def __init__(self, formula, dop_condition, decision=None):
        self.formula = formula                          # выражение вида y'(x) = f(x, y)
        self.dop_condition = dop_condition              # выражение вида y(0) = 1
        self.decision = decision                        # аналитическое решение уравнения

        self.derivative = self.formula[self.formula.find('=') + 1:].strip()
        self.x_point = float(self.dop_condition[self.dop_condition.find('(') + 1: self.dop_condition.find(')')].strip())
        self.y_point = float(self.dop_condition[self.dop_condition.find('=') + 1:].strip())

    def __str__(self):
        return f'Формула     : {self.formula}\n' \
               f'Доп условие : {self.dop_condition}\n' \
               f'Производная : {self.derivative}\n' \
               f'Х           : {self.x_point} \n' \
               f'Y           : {self.y_point}\n' \
               f'D           : {self.decision} \n'

    def get_decision_points(self, a, b, n):
        def if_y(y_a, y_b, n):
            """ Вывести граф по формуле x = ... """
            formula = self.decision[self.decision.find('=') + 1:].strip()

            _points = []
            y = y_a
            h = (y_b - y_a) / n

            while y <= y_b:
                x = eval(formula, {"y": y, "e": e})
                _points.append((x, y))

                y += h

            return _points

        def if_x(x_a, x_b, n):
            """ Вывести граф по формуле y = ... """
            formula = self.decision[self.decision.find('=') + 1:].strip()
            _points = []
            x = x_a
            h = (x_b - x_a) / n

            while x <= x_b:
                y = eval(formula, {"x": x})
                _points.append((x, y))

                x += h

            return _points

        if self.decision.strip()[0] == 'x':
            return if_y(a, b, n)

        return if_x(a, b, n)


class Euler(ODY):
    def __init__(self, formula, dop_condition, decision=None):
        super().__init__(formula, dop_condition, decision)

        self.explicit_solution = None                   # Явное решение
        self.implicit_solution = None                   # Неявное решение

    def find_explicit_euler_decision(self, a, b, n):
        # a, b -- начало и конец отрезка
        # n -- количество узлов на сетке

        def left():
            x = self.x_point
            y = self.y_point

            while x >= a:
                self.explicit_solution.append((x, y))

                d = eval(self.derivative, {"x": x, "y": y})
                x = x - h
                y = y - h * d

        def right():
            x = self.x_point
            y = self.y_point

            while x <= b:
                self.explicit_solution.append((x, y))
                d = eval(self.derivative, {"x": x, "y": y})
                x = x + h
                y = y + h * d

        self.explicit_solution = []

        h = (b - a) / n

        if a <= self.x_point <= b:
            left()
            right()
        elif a <= self.x_point and b <= self.x_point:
            a, b = min(a, b), max(a, b)
            left()
            self.explicit_solution = [i for i in self.explicit_solution if i[0] <= b]
        elif a >= self.x_point and b >= self.x_point:
            a, b = min(a, b), max(a, b)
            right()
            self.explicit_solution = [i for i in self.explicit_solution if i[0] >= a]

        return self.explicit_solution


ody_1 = Euler("y'(x) = 1 / (y * y + x)", "y(1) = 0", "x = 3 * e ** y - y ** 2 - 2 * y - 2")

# points = ody_1.get_decision_points(-2, 1, 20)
# explicit_euler_points = ody_1.find_explicit_euler_decision(-5, 5, 20)

points = ody_1.get_decision_points(-5, 5, 50)
explicit_euler_points = ody_1.find_explicit_euler_decision(-5, 5, 20)

Graph.show_list(points, label='analytical')
Graph.show_list(explicit_euler_points, 'm', label='explicit Euler', show_points=True)

Graph.draw_all()