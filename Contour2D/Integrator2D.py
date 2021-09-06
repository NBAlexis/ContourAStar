import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour2D.Integrand2D import Integrand2D
from UsefulFunctions.GaussianPatterson import GaussPatterson
from UsefulFunctions.NestedQuadrature import NestedQuadrature
from UsefulFunctions.SparseGridGenerator import SparseGrid


class Integrators2D:
    def Integrate(self, func: Integrand2D,
                  fromX: complex, toX: complex, fromY: complex, toY: complex) -> [bool, complex]:
        pass

    def PartialIntegrateY(self, func: Integrand2D, x: complex, fromY: complex, toY: complex) -> [bool, complex]:
        pass

    def GetLeftEdgeX(self) -> float:
        pass


class SparseGridIntegrator(Integrators2D):

    def __init__(self, nestedQuadrature: NestedQuadrature = GaussPatterson(), maxOrder: int = 15, epsilon=1.0e-8,
                 epsilon1d=1.0e-8, checkContinuity: float = -1, logLevel: LogLevel = LogLevel.General):
        self.epsilon = epsilon
        self.epsilon1d = epsilon1d
        self.logLevel = logLevel
        sparseGrid = SparseGrid(nestedQuadrature)
        self.maxOrder = maxOrder if nestedQuadrature.MaxOrder() < 0 else min(maxOrder, nestedQuadrature.MaxOrder())
        [pts, startIndex, endIndex, weights] = sparseGrid.ConstructPointListAndWeightList(self.maxOrder)
        self.points = pts
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.weights = weights
        [pts2, startIndex2, endIndex2, weights2] = nestedQuadrature.ConstructPointListYLeftMost(self.maxOrder)
        self.pointsY = pts2
        self.startIndexY = startIndex2
        self.endIndexY = endIndex2
        self.weightsY = weights2
        self.checkContinuity = checkContinuity
        self.nestedQuadrature = nestedQuadrature

    def Integrate(self, func: Integrand2D,
                  fromX: complex, toX: complex, fromY: complex, toY: complex) -> [bool, complex]:
        if self.logLevel >= LogLevel.Verbose:
            print("integrate: x: ", fromX, " to ", toX, " y: ", fromY, " to ", toY)
        resOld: complex = 0
        strideX: complex = 0.5 * (toX - fromX)
        strideY: complex = 0.5 * (toY - fromY)
        delta = 0
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            for pointIndex in range(self.startIndex[i], self.endIndex[i]):
                x = self.points[pointIndex].x
                y = self.points[pointIndex].y
                v = func.Evaluate(strideX * (x + 1) + fromX, strideY * (y + 1) + fromY)
                if cmath.isnan(v) or cmath.isinf(v):
                    print("Sparse Grid failed because of nan at {}, {}".format(strideX * (x + 1) + fromX,
                                                                               strideY * (y + 1) + fromY))
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
                    return [True, strideX * strideY * resNew]
            resOld = resNew
        if self.logLevel >= LogLevel.General:
            print("Sparse Grid failed due to max iteration reached\n, x:{} to {} y:{} to {} last value is "
                    .format(fromX, toX, fromY, toY),
                    strideX * strideY * resOld,
                " delta = {}/{}".format(delta, self.epsilon))
        return [False, strideX * strideY * resOld]

    def PartialIntegrateY(self, func: Integrand2D, x: complex, fromY: complex, toY: complex) -> [bool, complex]:
        if self.logLevel > LogLevel.Verbose:
            print("integrate1D: y: ", fromY, " to ", toY)
        resOld: complex = 0
        strideY: complex = 0.5 * (toY - fromY)
        delta = 0
        # print(self.startIndexY)
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            for pointIndex in range(self.startIndexY[i], self.endIndexY[i]):
                y = self.pointsY[pointIndex].y
                v = func.Evaluate(x, strideY * (y + 1) + fromY)
                if cmath.isnan(v) or cmath.isinf(v):
                    print("Sparse Grid failed because of nan at  {}, {}, ({})".format(x, strideY * (y + 1) + fromY, y))
                    return [False, cmath.nan]
                if self.checkContinuity > 0:
                    dyReal = abs(func.Evaluate(x, strideY * (y + 1) + fromY + 1.0e-8) - func.Evaluate(x, strideY * (
                            y + 1) + fromY - 1.0e-8))
                    if dyReal > self.checkContinuity:
                        print("Sparse Grid failed because of discontinuity(R) at  {}, {}, ({})"
                              .format(x, strideY * (y + 1) + fromY, y))
                        return [False, cmath.nan]
                    dyImag = abs(func.Evaluate(x, strideY * (y + 1) + fromY + 1.0e-8j) - func.Evaluate(x, strideY * (
                            y + 1) + fromY - 1.0e-8j))
                    if dyImag > self.checkContinuity:
                        print("Sparse Grid failed because of discontinuity(I) at  {}, {}, ({})"
                              .format(x, strideY * (y + 1) + fromY, y))
                self.pointsY[pointIndex].v = v
            # ======= time weights ===============
            resNew: complex = 0
            for pointIndex in range(0, self.endIndexY[i]):
                resNew = resNew + self.pointsY[pointIndex].v * self.weightsY[i][pointIndex]
            if i > 0:
                checkDelta = abs(resNew)
                checkDelta = 1.0e-6 if checkDelta < 1.0e-6 else checkDelta
                delta = abs((resNew - resOld) / checkDelta)
                if self.logLevel > LogLevel.Verbose:
                    print("x = {}, Y from {} to {} = {}, delta = {}/{}"
                          .format(x, fromY, toY, strideY * resNew, delta, self.epsilon1d))
                if delta < self.epsilon1d:
                    if self.logLevel >= LogLevel.Verbose:
                        print("x = {}, Y from {} to {} = {}".format(x, fromY, toY, strideY * resNew))
                    return [True, strideY * resNew]
            resOld = resNew
        print(
            "Sparse Grid PartialIntegrateYEdge failed due to max iteration reached\n x={}, y: {} to {}, last value is ".format(
                x, fromY, toY), strideY * resOld, " last delta: {}/{}".format(delta, self.epsilon1d))
        return [False, strideY * resOld]

    def GetLeftEdgeX(self) -> float:
        [v, _, _] = self.nestedQuadrature.GetLeftMostPoint()
        return v
