import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour3D.Integrand3D import Integrand3D
from UsefulFunctions.GaussianPatterson import GaussPatterson
from UsefulFunctions.NestedQuadrature import NestedQuadrature
from UsefulFunctions.SparseGridGenerator import SparseGrid
from UsefulFunctions.SparseGridGenerator3D import SparseGrid3D


class Integrators3D:
    def Integrate(self,
                  func: Integrand3D,
                  fromX: complex,
                  toX: complex,
                  fromY: complex,
                  toY: complex,
                  fromZ: complex,
                  toZ: complex) -> [bool, complex]:
        pass

    def IntegrateYZ(self, func: Integrand3D,
                    x: complex,
                    fromY: complex,
                    toY: complex,
                    fromZ: complex,
                    toZ: complex) -> [bool, complex]:
        pass

    def IntegrateZ(self, func: Integrand3D,
                   x: complex,
                   y: complex,
                   fromZ: complex,
                   toZ: complex) -> [bool, complex]:
        pass

    def GetLeftEdgeX(self) -> float:
        pass


class SparseGridIntegrator3D(Integrators3D):

    def __init__(self, nestedQuadrature: NestedQuadrature = GaussPatterson(),
                 maxOrder: int = 15, epsilon=1.0e-8, epsilon2d=1.0e-8, epsilon1d=1.0e-8,
                 checkContinuity: float = -1, logLevel: LogLevel = LogLevel.General):
        self.epsilon = epsilon
        self.epsilon2d = epsilon2d
        self.epsilon1d = epsilon1d
        self.logLevel = logLevel
        self.maxOrder = maxOrder if nestedQuadrature.MaxOrder() < 0 else min(maxOrder, nestedQuadrature.MaxOrder())
        self.checkContinuity = checkContinuity
        self.nestedQuadrature = nestedQuadrature
        sparseGrid = SparseGrid3D(nestedQuadrature)
        [pts, startIndex, endIndex, weights] = sparseGrid.ConstructPointListAndWeightList(self.maxOrder)
        self.points = pts
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.weights = weights
        sparseGrid2d = SparseGrid(nestedQuadrature)
        [pts2, startIndex2, endIndex2, weights2] = sparseGrid2d.ConstructPointListAndWeightList(self.maxOrder)
        self.points2d = pts2
        self.startIndex2d = startIndex2
        self.endIndex2d = endIndex2
        self.weights2d = weights2
        [pts3, startIndex3, endIndex3, weights3] = nestedQuadrature.ConstructPointListYLeftMost(self.maxOrder)
        self.points1d = pts3
        self.startIndex1d = startIndex3
        self.endIndex1d = endIndex3
        self.weights1d = weights3

    def Integrate(self, func: Integrand3D,
                  fromX: complex, toX: complex,
                  fromY: complex, toY: complex,
                  fromZ: complex, toZ: complex) -> [bool, complex]:
        if self.logLevel >= LogLevel.Verbose:
            print("integrate3d: x: ", fromX, " to ", toX, " y: ", fromY, " to ", toY, " z: ", fromZ, " to ", toZ)
        resOld: complex = 0
        strideX: complex = 0.5 * (toX - fromX)
        strideY: complex = 0.5 * (toY - fromY)
        strideZ: complex = 0.5 * (toZ - fromZ)
        delta = 0
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            for pointIndex in range(self.startIndex[i], self.endIndex[i]):
                x = self.points[pointIndex].x
                y = self.points[pointIndex].y
                z = self.points[pointIndex].z
                v = func.Evaluate(strideX * (x + 1) + fromX, strideY * (y + 1) + fromY, strideZ * (z + 1) + fromZ)
                if cmath.isnan(v) or cmath.isinf(v):
                    print("Sparse Grid3d failed because of nan at {}, {}, {}"
                          .format(strideX * (x + 1) + fromX,
                                  strideY * (y + 1) + fromY,
                                  strideZ * (z + 1) + fromZ))
                    return [False, cmath.nan]
                self.points[pointIndex].v = v
            # ======= time weights ===============
            resNew: complex = 0
            for pointIndex in range(0, self.endIndex[i]):
                resNew = resNew + self.points[pointIndex].v * self.weights[i][pointIndex]
            if i > 0:
                checkDelta = abs(resNew)
                checkDelta = 1.0e-6 if checkDelta < 1.0e-6 else checkDelta
                delta = abs((resNew - resOld) / checkDelta)
                if self.logLevel >= LogLevel.Verbose:
                    print("delta = ", delta)
                if delta < self.epsilon:
                    if self.logLevel >= LogLevel.Verbose:
                        print("Sparse Grid 3d done, value is ", strideX * strideY * strideZ * resNew)
                    return [True, strideX * strideY * strideZ * resNew]
            resOld = resNew
        if self.logLevel >= LogLevel.General:
            print("Sparse Grid 3d failed due to max iteration reached\n, x:{} to {} y:{} to {} z:{} to {}  last value is "
                  .format(fromX, toX, fromY, toY, fromZ, toZ),
                  strideX * strideY * strideZ * resOld,
                  " delta = {}/{}".format(delta, self.epsilon))
        return [False, strideX * strideY * strideZ * resOld]

    def IntegrateYZ(self, func: Integrand3D,
                    x: complex,
                    fromY: complex, toY: complex,
                    fromZ: complex, toZ: complex) -> [bool, complex]:
        if self.logLevel >= LogLevel.Verbose:
            print("integrate 2d: x: ", x, " y: ", fromY, " to ", toY, " z: ", fromZ, " to ", toZ)
        resOld: complex = 0
        strideY: complex = 0.5 * (toY - fromY)
        strideZ: complex = 0.5 * (toZ - fromZ)
        delta = 0
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            for pointIndex in range(self.startIndex2d[i], self.endIndex2d[i]):
                y = self.points2d[pointIndex].x
                z = self.points2d[pointIndex].y
                v = func.Evaluate(x, strideY * (y + 1) + fromY, strideZ * (z + 1) + fromZ)
                if cmath.isnan(v) or cmath.isinf(v):
                    print("Sparse Grid 2d failed because of nan at {}, {}, {}"
                          .format(x,
                                  strideY * (y + 1) + fromY,
                                  strideZ * (z + 1) + fromZ))
                    return [False, cmath.nan]
                self.points2d[pointIndex].v = v
            # ======= time weights ===============
            resNew: complex = 0
            for pointIndex in range(0, self.endIndex2d[i]):
                resNew = resNew + self.points2d[pointIndex].v * self.weights2d[i][pointIndex]
            if i > 0:
                checkDelta = abs(resNew)
                checkDelta = 1.0e-6 if checkDelta < 1.0e-6 else checkDelta
                delta = abs((resNew - resOld) / checkDelta)
                if self.logLevel >= LogLevel.Verbose:
                    print("delta = ", delta)
                if delta < self.epsilon2d:
                    return [True, strideY * strideZ * resNew]
            resOld = resNew
        if self.logLevel >= LogLevel.General:
            print("Sparse Grid 2d failed due to max iteration reached\n, x:{} y:{} to {} z:{} to {}  last value is "
                  .format(x, fromY, toY, fromZ, toZ),
                  strideY * strideZ * resOld,
                  " delta = {}/{}".format(delta, self.epsilon2d))
        return [False, strideY * strideZ * resOld]

    def IntegrateZ(self, func: Integrand3D,
                   x: complex, y: complex,
                   fromZ: complex, toZ: complex) -> [bool, complex]:
        if self.logLevel >= LogLevel.Verbose:
            print("integrate 1d: x: ", x, " y: ", y, " z: ", fromZ, " to ", toZ)
        resOld: complex = 0
        strideZ: complex = 0.5 * (toZ - fromZ)
        delta = 0
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            for pointIndex in range(self.startIndex1d[i], self.endIndex1d[i]):
                z = self.points1d[pointIndex].y
                v = func.Evaluate(x, y, strideZ * (z + 1) + fromZ)
                if cmath.isnan(v) or cmath.isinf(v):
                    print("Sparse Grid 1d failed because of nan at {}, {}, {}"
                          .format(x, y,
                                  strideZ * (z + 1) + fromZ))
                    return [False, cmath.nan]
                self.points1d[pointIndex].v = v
            # ======= time weights ===============
            resNew: complex = 0
            for pointIndex in range(0, self.endIndex1d[i]):
                resNew = resNew + self.points1d[pointIndex].v * self.weights1d[i][pointIndex]
            if i > 0:
                checkDelta = abs(resNew)
                checkDelta = 1.0e-6 if checkDelta < 1.0e-6 else checkDelta
                delta = abs((resNew - resOld) / checkDelta)
                if self.logLevel >= LogLevel.Verbose:
                    print("delta = ", delta)
                if delta < self.epsilon1d:
                    return [True, strideZ * resNew]
            resOld = resNew
        if self.logLevel >= LogLevel.General:
            print("Sparse Grid 1d failed due to max iteration reached\n, x:{} y:{} z:{} to {}  last value is "
                  .format(x, y, fromZ, toZ),
                  strideZ * resOld,
                  " delta = {}/{}".format(delta, self.epsilon1d))
        return [False, strideZ * resOld]

    def GetLeftEdgeX(self) -> float:
        [v, _, _] = self.nestedQuadrature.GetLeftMostPoint()
        return v
