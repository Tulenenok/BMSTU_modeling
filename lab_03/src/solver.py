"""
    Основная часть программы (итерационный метод)
"""

from dataclasses import dataclass
from typing import Callable
from typing import List

from runthrough import (
    run_through_method,
    Matrix3Diagonal,
    InitialValues,
)
from tools import (
    TableFunc,
    trapezoidal
)


@dataclass
class Input:
    lambda_f: Callable[[float], float]
    k_f: Callable[[float], float]
    n: float = 0
    r0: float = 0
    r_max: float = 0
    steps: int = 0
    t0: float = 0
    sigma: float = 0
    f_0: float = 0
    alpha: float = 0
    eps1: float = 1e-4
    eps2: float = 1e-4


@dataclass
class Output:
    radius: List[float]
    temperature: List[float]
    f1: float
    f2: float
    iter: int


def solve(params: Input) -> Output:
    """
        Используя метод прогонки на каждой итерации уточняем коэффициенты и
        приближаемся к правильному ответу
    """

    step = (params.r_max - params.r0) / params.steps
    ir_max = 1 / params.r_max
    origin_grid = [params.r0 + i * step for i in range(params.steps + 1)]
    grid = [i * ir_max for i in origin_grid]
    res = [params.t0 for _ in origin_grid]

    irsqrh = 1 / (params.r_max * params.r_max * step)
    ir = 1 / params.r_max
    basep = 4 * params.n * params.n * params.sigma
    basef = basep * params.t0 ** 4
    at0ir = params.alpha * params.t0 / params.r_max

    matrix = Matrix3Diagonal(params.steps - 1,
                             [0] * (params.steps - 1),
                             [0] * (params.steps - 1),
                             [0] * (params.steps - 1),
                             [0] * (params.steps - 1))
    start = InitialValues()
    end = InitialValues()
    f1 = 0
    f2 = 0

    run = True
    iter = 0

    while (run):
        for i in range(matrix.size):
            zm = (grid[i] + grid[i + 1]) / 2
            zp = (grid[i + 2] + grid[i + 1]) / 2
            v = (zp * zp - zm * zm) / 2
            kappam = (params.lambda_f(res[i])
                      + params.lambda_f(res[i + 1])) / 2
            kappap = (params.lambda_f(res[i + 2])
                      + params.lambda_f(res[i + 1])) / 2

            k = params.k_f(res[i + 1])
            p = basep * k * res[i + 1] ** 3
            f = basef * k

            matrix.A[i] = irsqrh * zm * kappam
            matrix.C[i] = irsqrh * zp * kappap
            matrix.B[i] = matrix.A[i] + matrix.C[i] + p * v
            matrix.D[i] = f * v

        zhl = (grid[0] + grid[1]) / 2
        zhh = (grid[params.steps] + grid[params.steps - 1]) / 2

        kappal = (params.lambda_f(res[0]) + params.lambda_f(res[1])) / 2
        kappah = (params.lambda_f(res[params.steps])
                  + params.lambda_f(res[params.steps - 1])) / 2

        k0 = params.k_f(res[0])
        k1 = params.k_f(res[1])
        kN0 = params.k_f(res[params.steps - 1])
        kN1 = params.k_f(res[params.steps])

        p0 = basep * k0 * res[0] ** 3
        p1 = basep * k1 * res[1] ** 3
        pN0 = basep * kN0 * res[params.steps - 1] ** 3
        pN1 = basep * kN1 * res[params.steps] ** 3

        pl = (p0 + p1) / 2
        ph = (pN0 + pN1) / 2

        f0 = basef * k0
        f1 = basef * k1
        fN0 = basef * kN0
        fN1 = basef * kN1

        fl = (f0 + f1) / 2
        fh = (fN0 + fN1) / 2

        start.m = (kappal * irsqrh + step * pl / 8) * zhl + step * p0 * grid[0] / 4
        start.k = (step * pl / 8 - kappal * irsqrh) * zhl
        start.p = grid[0] * params.f_0 * ir + step * (f0 * grid[0] + fl * zhl) / 4

        end.m = (-kappah * irsqrh + step * ph / 8) * zhh
        end.k = kappah * zhh * irsqrh + grid[params.steps] * params.alpha * ir \
                + step * (ph * zhh / 2 + pN1 * grid[params.steps]) / 4
        end.p = at0ir * grid[params.steps] + \
                + step * (fN1 * grid[params.steps] + fh * zhh) / 4

        res_new = run_through_method(matrix, start, end)
        c = 0

        for i in range(params.steps + 1):
            if abs((res[i] - res_new[i]) / res_new[i]) < params.eps1:
                c += 1

        if c == params.steps + 1:
            run = False

        if not run:
            f1 = params.r0 * params.f_0 \
                 - params.r_max * params.alpha * (res_new[-1] - params.t0)

            tmp = [params.k_f(res_new[i]) * (res_new[i] ** 4
                                             - params.t0 ** 4)
                   * origin_grid[i] for i in range(params.steps + 1)]

            f2 = 4 * params.n * params.n * params.sigma \
                 * trapezoidal(params.r0, params.r_max, params.steps,
                               TableFunc(origin_grid, tmp))

            if 1e-8 > abs(f1):
                if 1e-8 < abs(f2):
                    run = True
            elif abs((f1 - f2) / f1) > params.eps2:
                run = True

        res = res_new
        iter += 1

    return Output(origin_grid, res, f1, f2, iter)

