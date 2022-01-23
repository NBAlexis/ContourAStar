import cmath
import random

import numpy as np

from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import Integrators3D


class MCIntegrator3D(Integrators3D):

    def __init__(self, pointCount: int = 10000, epsilon1: float = 0.01, epsilon2: float = 0.2, smallProtect: float = 1.0e-15):
        self.pointCount = pointCount
        self.epsilon1 = epsilon1
        self.epsilon2 = epsilon2
        self.point1 = np.random.rand(pointCount, 3)
        self.point100 = np.random.rand(pointCount * 100, 3)
        self.point10000 = np.random.rand(pointCount * 10000, 3)
        self.smallProtect = smallProtect

    def Integrate(self, func: Integrand3D,
                  fromX: complex, toX: complex,
                  fromY: complex, toY: complex,
                  fromZ: complex, toZ: complex) -> [bool, complex]:
        pfx = lambda a: a * (toX - fromX) + fromX
        pfy = lambda a: a * (toY - fromY) + fromY
        pfz = lambda a: a * (toZ - fromZ) + fromZ
        vpfx = np.vectorize(pfx)
        vpfy = np.vectorize(pfy)
        vpfz = np.vectorize(pfz)
        p1x = vpfx(self.point1[:, 0])
        p1y = vpfy(self.point1[:, 1])
        p1z = vpfz(self.point1[:, 2])
        p100x = vpfx(self.point100[:, 0])
        p100y = vpfy(self.point100[:, 1])
        p100z = vpfz(self.point100[:, 2])
        p10000x = vpfx(self.point10000[:, 0])
        p10000y = vpfy(self.point10000[:, 1])
        p10000z = vpfz(self.point10000[:, 2])
        vpf = func.GetVPF()
        try:
            v1 = vpf(p1x, p1y, p1z)
            v100 = vpf(p100x, p100y, p100z)
            v10000 = vpf(p10000x, p10000y, p10000z)
        except (ValueError, ZeroDivisionError):
            print("result has nan or inf")
            return [False, cmath.nan]
        mean1 = np.mean(v1)
        mean100 = np.mean(v100)
        mean10000 = np.mean(v10000)
        res = (mean1 * self.pointCount + mean100 * self.pointCount * 100 + mean10000 * self.pointCount * 10000) / (self.pointCount * 10101)
        if cmath.isnan(res) or cmath.isinf(res):
            print("result has nan or inf")
            return [False, cmath.nan]
        delta1 = abs((mean10000 - mean100) / mean10000)
        if delta1 > self.epsilon1:
            print("poor convergence: res = {}, delta = {} / {}".format(res, delta1, self.epsilon1))
            return [False, res]
        delta2 = abs((mean10000 - mean100) / (mean100 - mean1))
        if delta2 > self.epsilon2:
            print("poor convergence: res = {}, delta = {} / {} > epsilon: {}".format(res, abs(mean10000 - mean100), abs(mean100 - mean1), self.epsilon2))
            return [False, res]
        return [True, res]

    def IntegrateYZ(self, func: Integrand3D,
                    x: complex,
                    fromY: complex, toY: complex,
                    fromZ: complex, toZ: complex) -> [bool, complex]:
        pfy = lambda a: a * (toY - fromY) + fromY
        pfz = lambda a: a * (toZ - fromZ) + fromZ
        vpfy = np.vectorize(pfy)
        vpfz = np.vectorize(pfz)
        p1y = vpfy(self.point1[:, 1])
        p1z = vpfz(self.point1[:, 2])
        p100y = vpfy(self.point100[:, 1])
        p100z = vpfz(self.point100[:, 2])
        p10000y = vpfy(self.point10000[:, 1])
        p10000z = vpfz(self.point10000[:, 2])
        vpf = func.GetVPF()
        try:
            v1 = vpf(x, p1y, p1z)
            v100 = vpf(x, p100y, p100z)
            v10000 = vpf(x, p10000y, p10000z)
        except (ValueError, ZeroDivisionError):
            print("result has nan or inf")
            return [False, cmath.nan]
        mean1 = np.mean(v1)
        mean100 = np.mean(v100)
        mean10000 = np.mean(v10000)
        res = (mean1 * self.pointCount + mean100 * self.pointCount * 100 + mean10000 * self.pointCount * 10000) / (self.pointCount * 10101)
        if cmath.isnan(res) or cmath.isinf(res):
            print("result has nan or inf")
            return [False, cmath.nan]
        delta1 = abs((mean10000 - mean100) / mean10000)
        if delta1 > self.epsilon1:
            print("poor convergence: res = {}, delta = {} / {}".format(res, delta1, self.epsilon1))
            return [False, res]
        delta2 = abs((mean10000 - mean100) / (mean100 - mean1))
        if delta2 > self.epsilon2:
            print("poor convergence: res = {}, delta = {} / {} > epsilon: {}".format(res, abs(mean10000 - mean100), abs(mean100 - mean1), self.epsilon2))
            return [False, res]
        return [True, res]

    def IntegrateZ(self, func: Integrand3D,
                   x: complex, y: complex,
                   fromZ: complex, toZ: complex) -> [bool, complex]:
        pfz = lambda a: a * (toZ - fromZ) + fromZ
        vpfz = np.vectorize(pfz)
        p1z = vpfz(self.point1[:, 2])
        p100z = vpfz(self.point100[:, 2])
        p10000z = vpfz(self.point10000[:, 2])
        vpf = func.GetVPF()
        try:
            v1 = vpf(x, y, p1z)
            v100 = vpf(x, y, p100z)
            v10000 = vpf(x, y, p10000z)
        except (ValueError, ZeroDivisionError):
            print("result has nan or inf")
            return [False, cmath.nan]
        mean1 = np.mean(v1)
        mean100 = np.mean(v100)
        mean10000 = np.mean(v10000)
        res = (mean1 * self.pointCount + mean100 * self.pointCount * 100 + mean10000 * self.pointCount * 10000) / (self.pointCount * 10101)
        if cmath.isnan(res) or cmath.isinf(res):
            print("result has nan or inf")
            return [False, cmath.nan]
        delta1 = abs((mean10000 - mean100) / mean10000)
        if delta1 > self.epsilon1:
            print("poor convergence: res = {}, delta = {} / {}".format(res, delta1, self.epsilon1))
            return [False, res]
        delta2 = abs((mean10000 - mean100) / (mean100 - mean1))
        if delta2 > self.epsilon2:
            print("poor convergence: res = {}, delta = {} / {} > epsilon: {}".format(res, abs(mean10000 - mean100), abs(mean100 - mean1), self.epsilon2))
            return [False, res]
        return [True, res]

    def GetLeftEdgeX(self) -> float:
        return -1 + self.smallProtect
