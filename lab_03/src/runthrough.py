"""
    Реализация метода прогонки
"""

from dataclasses import dataclass
from typing import List


@dataclass
class Matrix3Diagonal:
    size: int
    A: List[float]
    B: List[float]
    C: List[float]
    D: List[float]


@dataclass
class InitialValues:
    m: float = 0
    k: float = 0
    p: float = 0


@dataclass
class RunThroughCoefficient:
    k: float = 0
    b: float = 0


def run_through_method(matrix, start: InitialValues, end: InitialValues):

    result = [0.0] * (matrix.size + 2)
    k = RunThroughCoefficient(-(start.k / start.m), start.p / start.m)
    run_coefficients: list[RunThroughCoefficient] = [k]

    for i in range(matrix.size):
        k = RunThroughCoefficient(
            matrix.C[i] / (matrix.B[i] - matrix.A[i] * k.k),
            (matrix.A[i] * k.b + matrix.D[i]) / (matrix.B[i] - matrix.A[i] * k.k)
        )

        run_coefficients.append(k)

    result[-1] = (end.p - end.m * k.b) / (end.k + end.m * k.k)

    for i in range(matrix.size, -1, -1):
        result[i] = run_coefficients[i].k * result[i + 1] \
                    + run_coefficients[i].b

    return result
