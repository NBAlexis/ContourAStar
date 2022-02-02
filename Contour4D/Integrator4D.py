import cmath
import math

from Contour1D.CommonDefinitions import LogLevel
from Contour4D.Integrand4D import Integrand4D
from SparseGridIntegrators.GaussianPatterson import GaussPatterson
from SparseGridIntegrators.NestedQuadrature import NestedQuadrature
from SparseGridIntegrators.SparseGridGenerator4D import SparseGrid4D


class Integrators4D:
    def Integrate(self,
                  func: Integrand4D,
                  fromX: complex,
                  toX: complex,
                  fromY: complex,
                  toY: complex,
                  fromZ: complex,
                  toZ: complex,
                  fromW: complex,
                  toW: complex) -> [bool, complex]:
        pass

    def Finish(self):
        pass

    @staticmethod
    def __FillLine(func: Integrand4D,
                   xf: complex,
                   xt: complex,
                   yf: complex,
                   yt: complex,
                   zf: complex,
                   zt: complex,
                   wf: complex,
                   wt: complex,
                   step: list,
                   interval: int,
                   d: int,
                   argLst: list) -> float:
        if 0 == d:
            sampLst = [func.EvaluateDenominator(xf + i * step[0],
                                                yf if 0 == argLst[0] else yt,
                                                zf if 0 == argLst[1] else zt,
                                                wf if 0 == argLst[2] else wt)
                       for i in range(0, interval + 1)]
            ret: float = 0.0
            for i in range(0, interval):
                ret = ret + cmath.phase(sampLst[i + 1] / sampLst[i])
            return ret
        if 1 == d:
            sampLst = [func.EvaluateDenominator(xf if 0 == argLst[0] else xt,
                                                yf + i * step[1],
                                                zf if 0 == argLst[1] else zt,
                                                wf if 0 == argLst[2] else wt)
                       for i in range(0, interval + 1)]
            ret: float = 0.0
            for i in range(0, interval):
                ret = ret + cmath.phase(sampLst[i + 1] / sampLst[i])
            return ret
        if 2 == d:
            sampLst = [func.EvaluateDenominator(xf if 0 == argLst[0] else xt,
                                                yf if 0 == argLst[1] else yt,
                                                zf + i * step[2],
                                                wf if 0 == argLst[2] else wt)
                       for i in range(0, interval + 1)]
            ret: float = 0.0
            for i in range(0, interval):
                ret = ret + cmath.phase(sampLst[i + 1] / sampLst[i])
            return ret
        sampLst = [func.EvaluateDenominator(xf if 0 == argLst[0] else xt,
                                            yf if 0 == argLst[1] else yt,
                                            zf if 0 == argLst[2] else zt,
                                            wf + i * step[3])
                   for i in range(0, interval + 1)]
        ret: float = 0.0
        for i in range(0, interval):
            ret = ret + cmath.phase(sampLst[i + 1] / sampLst[i])
        return ret

    @staticmethod
    def CheckPolesFast(func: Integrand4D,
                       fromX: complex,
                       toX: complex,
                       fromY: complex,
                       toY: complex,
                       fromZ: complex,
                       toZ: complex,
                       fromW: complex,
                       toW: complex,
                       interval: int = 50,
                       epsilon: float = 1.0e-9) -> bool:
        """
        :return:
        Whether has a pole

        To avoid treat zeros as the pole, please specify the denominator.
        """
        if not func.HasDenominator():
            return False
        # rescale X
        sepX = (toX - fromX)
        xf = fromX + sepX * epsilon
        xt = toX - sepX * epsilon
        xStep = (xt - xf) / interval
        # rescale Y
        sepY = (toY - fromY)
        yf = fromY + sepY * epsilon
        yt = toY - sepY * epsilon
        yStep = (yt - yf) / interval
        # rescale Z
        sepZ = (toZ - fromZ)
        zf = fromZ + sepZ * epsilon
        zt = toZ - sepZ * epsilon
        zStep = (zt - zf) / interval
        # rescale W
        sepW = (toW - fromW)
        wf = fromW + sepW * epsilon
        wt = toW - sepW * epsilon
        wStep = (wt - wf) / interval
        stepAll = [xStep, yStep, zStep, wStep]
        listArgs = [[[[math.nan for _ in range(0, 4)] for _ in range(0, 2)] for _ in range(0, 2)] for _ in range(0, 2)]
        for d1 in range(0, 4):
            for d2 in range(d1 + 1, 4):
                for p1 in range(0, 2):
                    for p2 in range(0, 2):
                        # This is d1-d2 face at (p1, p2)
                        # line 1: (d1, p1, p2, 0)
                        # line 2: (d2, p1, p2, 0)
                        lineArgs1 = []
                        lineArgs2 = []
                        lineArgs3 = []
                        lineArgs4 = []
                        p1putA = False
                        p1putB = False
                        for k in range(0, 4):
                            if k != d1:
                                if k == d2:
                                    lineArgs1.append(0)
                                    lineArgs3.append(1)
                                else:
                                    if not p1putA:
                                        lineArgs1.append(p1)
                                        lineArgs3.append(p1)
                                        p1putA = True
                                    else:
                                        lineArgs1.append(p2)
                                        lineArgs3.append(p2)
                            if k != d2:
                                if k == d1:
                                    lineArgs2.append(1)
                                    lineArgs4.append(0)
                                else:
                                    if not p1putB:
                                        lineArgs2.append(p1)
                                        lineArgs4.append(p1)
                                        p1putB = True
                                    else:
                                        lineArgs2.append(p2)
                                        lineArgs4.append(p2)
                        # argLine1: float = 0.0
                        # argLine2: float = 0.0
                        # argLine3: float = 0.0
                        # argLine4: float = 0.0
                        if math.isnan(listArgs[d1][lineArgs1[0]][lineArgs1[1]][lineArgs1[2]]):
                            argLine1 = Integrators4D.__FillLine(func, xf, xt, yf, yt, zf, zt, wf, wt,
                                                                stepAll, interval, d1, lineArgs1)
                            listArgs[d1][lineArgs1[0]][lineArgs1[1]][lineArgs1[2]] = argLine1
                        else:
                            argLine1 = listArgs[d1][lineArgs1[0]][lineArgs1[1]][lineArgs1[2]]
                        if math.isnan(listArgs[d2][lineArgs2[0]][lineArgs2[1]][lineArgs2[2]]):
                            argLine2 = Integrators4D.__FillLine(func, xf, xt, yf, yt, zf, zt, wf, wt,
                                                                stepAll, interval, d2, lineArgs2)
                            listArgs[d2][lineArgs2[0]][lineArgs2[1]][lineArgs2[2]] = argLine2
                        else:
                            argLine2 = listArgs[d2][lineArgs2[0]][lineArgs2[1]][lineArgs2[2]]
                        if math.isnan(listArgs[d1][lineArgs3[0]][lineArgs3[1]][lineArgs3[2]]):
                            argLine3 = Integrators4D.__FillLine(func, xf, xt, yf, yt, zf, zt, wf, wt,
                                                                stepAll, interval, d1, lineArgs3)
                            listArgs[d1][lineArgs3[0]][lineArgs3[1]][lineArgs3[2]] = argLine3
                        else:
                            argLine3 = listArgs[d1][lineArgs3[0]][lineArgs3[1]][lineArgs3[2]]
                        if math.isnan(listArgs[d2][lineArgs4[0]][lineArgs4[1]][lineArgs4[2]]):
                            argLine4 = Integrators4D.__FillLine(func, xf, xt, yf, yt, zf, zt, wf, wt,
                                                                stepAll, interval, d2, lineArgs4)
                            listArgs[d2][lineArgs4[0]][lineArgs4[1]][lineArgs4[2]] = argLine4
                        else:
                            argLine4 = listArgs[d2][lineArgs4[0]][lineArgs4[1]][lineArgs4[2]]
                        faceArg = abs(argLine1 + argLine2 - argLine3 - argLine4)
                        if faceArg > 0.5:
                            return True
        return False


class SparseGridIntegrator4D(Integrators4D):

    def __init__(self, nestedQuadrature: NestedQuadrature = GaussPatterson(),
                 maxOrder: int = 15, epsilon=1.0e-8,
                 fastCheckPole: int = 20, logLevel: LogLevel = LogLevel.General):
        self.epsilon = epsilon
        self.logLevel = logLevel
        self.maxOrder = maxOrder if nestedQuadrature.MaxOrder() < 0 else min(maxOrder, nestedQuadrature.MaxOrder())
        self.nestedQuadrature = nestedQuadrature
        sparseGrid = SparseGrid4D(nestedQuadrature)
        [pts, startIndex, endIndex, weights] = sparseGrid.ConstructPointListAndWeightList(self.maxOrder)
        self.points = pts
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.weights = weights
        self.fastCheckPole = fastCheckPole

    def Integrate(self, func: Integrand4D,
                  fromX: complex, toX: complex,
                  fromY: complex, toY: complex,
                  fromZ: complex, toZ: complex,
                  fromW: complex, toW: complex) -> [bool, complex]:
        if self.logLevel >= LogLevel.Verbose:
            print("integrate4d: x: ", fromX, " to ", toX, " y: ", fromY, " to ", toY, " z: ", fromZ, " to ", toZ,
                  " w: ", fromW, " to ", toW)
        if self.fastCheckPole > 0:
            if self.CheckPolesFast(func, fromX, toX, fromY, toY, fromZ, toZ, fromW, toW, self.fastCheckPole):
                return [False, cmath.nan]
        resOld: complex = 0
        strideX: complex = 0.5 * (toX - fromX)
        strideY: complex = 0.5 * (toY - fromY)
        strideZ: complex = 0.5 * (toZ - fromZ)
        strideW: complex = 0.5 * (toW - fromW)
        delta = 0
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            for pointIndex in range(self.startIndex[i], self.endIndex[i]):
                x = self.points[pointIndex].x
                y = self.points[pointIndex].y
                z = self.points[pointIndex].z
                w = self.points[pointIndex].w
                v = func.Evaluate(strideX * (x + 1) + fromX, strideY * (y + 1) + fromY, strideZ * (z + 1) + fromZ,
                                  strideW * (w + 1) + fromW)
                if cmath.isnan(v) or cmath.isinf(v):
                    print("Sparse Grid3d failed because of nan at {}, {}, {}, {}"
                          .format(strideX * (x + 1) + fromX,
                                  strideY * (y + 1) + fromY,
                                  strideZ * (z + 1) + fromZ,
                                  strideW * (w + 1) + fromW))
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
                        print("Sparse Grid 4d done, value is ", strideX * strideY * strideZ * strideW * resNew)
                    return [True, strideX * strideY * strideZ * strideW * resNew]
            resOld = resNew
        if self.logLevel >= LogLevel.General:
            print(
                "Sparse Grid 3d failed due to max iteration reached\n, x:{} to {} y:{} to {} z:{} to {} w:{} to {}  last value is "
                    .format(fromX, toX, fromY, toY, fromZ, toZ, fromW, toW),
                    strideX * strideY * strideZ * strideW * resOld,
                " delta = {}/{}".format(delta, self.epsilon))
        return [False, strideX * strideY * strideZ * strideW * resOld]

    def Finish(self):
        return
