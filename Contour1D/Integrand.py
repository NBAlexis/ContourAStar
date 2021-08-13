import cmath
from enum import IntEnum


class IntegrandType(IntEnum):
    ZeroOne = 0
    AB = 1
    ZeroInfinite = 2
    AInfinite = 3
    InfiniteZero = 4
    InfiniteA = 5
    InfiniteInfinite = 6


class Integrand:

    def __init__(self, func, left: complex, right: complex, case: IntegrandType = IntegrandType.AB, epsilon=1.0e-6,
                 maxIterate: int = 10, iterateStep=1e-15, tolerance=0.1):
        self.func = func
        self.case: IntegrandType = case
        if cmath.isinf(left) and cmath.isinf(right):
            self.case = IntegrandType.InfiniteInfinite
        elif cmath.isinf(left):
            self.case = IntegrandType.InfiniteA
        elif cmath.isinf(right):
            self.case = IntegrandType.AInfinite
        if 0 == left and 1 == right:
            self.case = IntegrandType.ZeroOne
        self.left: complex = left
        self.right: complex = right
        self.sep: complex = 0
        if case == IntegrandType.AB:
            self.sep = right - left
        self._tolerance = tolerance
        self._epsilon = epsilon
        self._maxIterate = maxIterate
        self._iterateStep = iterateStep
        print(self.case)

    def EvaluateRight(self, x: complex, r: complex) -> complex:
        try:
            return self._Evaluate(x)
        except (ValueError, ZeroDivisionError):
            start = 1
            sep = r - x
            lastV = 0
            try:
                delta = 0
                for i in range(0, self._maxIterate):
                    start = start * self._iterateStep
                    retV = self._Evaluate(x + start * sep)
                    if 0 == i:
                        lastV = retV
                    else:
                        delta = abs(lastV - retV)
                        if delta < self._epsilon:
                            return retV
                        lastV = retV
                if delta < self._tolerance:
                    return lastV
                return cmath.nan
            except (ValueError, ZeroDivisionError):
                return cmath.nan

    def EvaluateLeftRight(self, x: complex, left: complex, right: complex) -> complex:
        try:
            return self._Evaluate(x)
        except (ValueError, ZeroDivisionError):
            start = 1
            sepR = right - x
            sepL = left - x
            lastVR = 0
            lastVL = 0
            try:
                delta = 0
                for i in range(0, self._maxIterate):
                    start = start * self._iterateStep
                    retVR = self._Evaluate(x + start * sepR)
                    retVL = self._Evaluate(x + start * sepL)
                    if 0 == i:
                        lastVR = retVR
                        lastVL = retVL
                    else:
                        delta = abs(lastVR - retVR) + abs(lastVL - retVL) + abs(retVL - retVR)
                        if delta < self._epsilon:
                            return 0.5 * (retVL + retVR)
                        lastVR = retVR
                        lastVL = retVL
                if delta < self._tolerance:
                    return 0.5 * (lastVR + lastVL)
                return cmath.nan
            except (ValueError, ZeroDivisionError):
                return cmath.nan

    def _Evaluate(self, x):
        if self.case == IntegrandType.AB:
            """
            Integrate[f[x], {x, a, b}] = (b-a) Integrate[f[(b - a)x + a], {x, 0, 1}] 
            """
            return self.sep * self.func(self.sep * x + self.left)
        elif self.case == IntegrandType.ZeroInfinite:
            """
            Integrate[f[x], {x, 0, inf}] = (1/x) Integrate[f[-Log[x]], {x, 0, 1}] 
            """
            return self.func(-cmath.log(x)) / x
        elif self.case == IntegrandType.InfiniteZero:
            """
            Integrate[f[x], {x, -inf, 0}] = (1/x) Integrate[f[Log[x]], {x, 0, 1}] 
            """
            return self.func(cmath.log(x)) / x
        elif self.case == IntegrandType.AInfinite:
            """
            Integrate[f[x], {x, a, inf}] = (1/x) Integrate[f[-Log[x] + a], {x, 0, 1}] 
            """
            return self.func(-cmath.log(x) + self.left) / x
        elif self.case == IntegrandType.InfiniteA:
            """
            Integrate[f[x], {x, -inf, a}] = (1/x) Integrate[f[Log[x] + a], {x, 0, 1}] 
            """
            return self.func(cmath.log(x) + self.right) / x
        elif self.case == IntegrandType.InfiniteInfinite:
            """
            Integrate[f[x], {x, -inf, inf}] = Integrate[ pi sec^2 [pi (x - 1/2)] f[ tan[(pi (x - 1/2))] ], {x, 0, 1} ] 
            """
            b: complex = cmath.pi * (x - 0.5)
            return cmath.pi * self.func(cmath.tan(b)) / (cmath.cos(b) ** 2)
        return self.func(x)
