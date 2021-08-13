import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour1D.Integrand import Integrand


class Integrators:
    def Integrate(self, func: Integrand, fromV: complex, toV: complex) -> [bool, complex]:
        pass


class Simpson(Integrators):

    def __init__(self, epsilon=1.0e-6, maxIterate: int = 15, logLevel: LogLevel = LogLevel.General,
                 Romberg: bool = True):
        self.epsilon = epsilon
        self.maxIterate = maxIterate
        self.logLevel = logLevel
        self.Romberg = Romberg

    def Integrate(self, func: Integrand, fromV: complex, toV: complex) -> [bool, complex]:
        stride: complex = toV - fromV
        startEnd: complex = func.EvaluateRight(fromV, toV) + func.EvaluateRight(toV, fromV)
        mid: complex = func.EvaluateLeftRight(fromV + 0.5 * stride, fromV, toV)
        resOld: complex = 0
        for seps in range(1, 1 + self.maxIterate):
            numberOfPoints: int = 1 << seps
            # 1->0.25, 0.75, 2->0.125, 0.375, 0.625, 0.875
            rStart = -0.5 / numberOfPoints
            rStride = 1 / numberOfPoints
            newMid: complex = 0
            for i in range(0, numberOfPoints):
                rStart = rStart + rStride
                newMid = newMid + func.EvaluateLeftRight(fromV + rStart * stride, fromV, toV)
            if 1 == seps:
                resOld = (startEnd + 2 * mid + 4 * newMid) * rStride * stride / 6
                mid = mid + newMid
            else:
                resNew = (startEnd + 2 * mid + 4 * newMid) * rStride * stride / 6
                if cmath.isnan(resNew):
                    if self.logLevel >= LogLevel.General:
                        print("Simpson failed due to nan")
                    return [False, cmath.nan]
                delta = resNew - resOld
                if self.logLevel >= LogLevel.Verbose:
                    print("Simpson step: {}, value: {}, delta: {}".format(seps, resNew, delta))
                if abs(delta) < self.epsilon:
                    if self.Romberg:
                        return [True, (16 * resNew - resOld) / 15]
                    return [True, resNew]
                resOld = resNew
                mid = mid + newMid
        print("Simpson failed due to max iteration reached, last value is ", resOld)
        return [False, resOld]
