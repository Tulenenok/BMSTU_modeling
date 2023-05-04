"""
    Вспомогательные алгоритмы, использованные при решении поставленной задачи
    (интерполяция, вычисление интеграла)
"""

import sys
from typing import (
    Callable,
    List,
)


def trapezoidal(s: float, e: float, n: int, f: Callable[[float], float]) -> float:
    """
        Решение интеграла методом трапеций
          - s -- нижний предел интегрирования
          - e -- верхний предел интегрирования
          - n -- кол-во узлов на сетке
          - f -- подынтегральная функция
        Основная идея: разобьем интервал интегрирования на n интервалов, на каждом интервале
        аппроксимируем функцию линейной функцией и найдем площадь получившейся трапеции
    """
    h = (e - s) / n
    x = s + h

    res = f(s) + f(e)

    for _ in range(1, n):
        res += 2 * f(x)
        x += h

    return res * h / 2


def binary_search(lst, value) -> int:
    l = 0
    r = len(lst) - 1

    while l < r:
        mid = (l + r) // 2

        if lst[mid] < value:
            l = mid + 1
        else:
            r = mid

    return l


def linear_interpolation(x: float, lx: List[float], ly: List[float]) -> float:
    """
        Линейная интерполяция
          - x  -- значение для оси абсцисс, для которого нам нужно предсказать значение y
          - lx -- список уже известных x
          - lx -- список уже известных y для x-ов
        Основная идея: есть две точки (x1, y1) и (x2, y2), мы хотим найти значение
        функции y в точке x, нахощейся между этими двумя точками. Для этого найдем уравнение
        прямой, проходящей через две точки и вычислим значение функции на отрезке между ними.
        y = y1 + (y2 - y1) * (x - x1) / (x2 - x1)
    """
    if len(lx) != len(ly):
        return 0

    if x < lx[0]:
        return ly[0]

    if x > lx[-1]:
        return ly[-1]

    i = binary_search(lx, x)

    if len(lx) <= i:
        return 0

    if lx[i] > x:
        i -= 1

    if sys.float_info.epsilon > abs(lx[i] - x):
        return ly[i]

    return ly[i] + (ly[i + 1] - ly[i]) * (x - lx[i]) / (lx[i + 1] - lx[i])


def newton_interpolation(x: float, lx: List[float], ly: List[float]) -> float:
    """
        Интерполяция Ньютона
          - x  -- значение для оси абсцисс, для которого нам нужно предсказать значение y
          - lx -- список уже известных x
          - lx -- список уже известных y для x-ов
        Основная идея: есть набор точек, мы хотим построить многочлен P(x), который
        проходит через все эти точки. Будем использовать раздельные разности.

        P(x) = f[x0] + (x - x0)f[x0, x1] + (x - x0)(x - x1)f[x0, x1, x2] + ... +
               (x - x0)(x - x1)...(x - xn-1)f[x0, x1, ..., xn], где

        f[xi] = yi
        f[xi, xi+1, ..., xj] = (f[xi+1, ..., xj] - f[xi, ..., xj-1]) / (xi - xj)
    """
    res = 0
    k = newton_get_coefficients(ly, lx)
    handle = 1

    for i in range(len(k)):
        res += handle * k[i]
        handle *= (x - lx[i])

    return res


def newton_get_coefficients(y: List[float], x: List[float]) -> List[float]:
    if len(x) != len(y):
        raise Exception

    k = y.copy()
    out = [k[0]]

    limit = len(x)
    i = 0
    while 1 != limit:

        j = 0
        while limit - 1 > j:
            k[j] = (k[j] - k[j + 1]) / (x[j] - x[j + i + 1])

            j += 1

        out.append(k[0])
        i += 1
        limit -= 1

    return out


class TableFunc:
    def __init__(self, x : List[float], y : List[float],
                 method : Callable[[float, List[float], List[float]], float]=linear_interpolation):
        self._x = x.copy()
        self._y = y.copy()
        self._method = method

    def __call__(self, x : float) -> float:
        return self._method(x, self._x, self._y)

