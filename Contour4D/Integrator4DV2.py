import cmath
import math

import numpy as np

from Contour1D.CommonDefinitions import LogLevel
from Contour1D.IntegrandV2 import IntegrandV2
from SparseGridIntegrators.SparseGridCache import CachedPointListAndWeightList


class Integrators4DV2:

    def Integrate(self,
                  func: IntegrandV2,
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
    def GenerateEdgeList(interval: int, epsilon: float = 1.0e-9) -> list:
        step = (1 - 2 * epsilon) / interval
        emptyList = []
        for d in range(0, 4):
            dlist = []
            for p1 in range(0, 2):
                p1list = []
                for p2 in range(0, 2):
                    p2list = []
                    for p3 in range(0, 2):
                        steplst = [epsilon + i * step for i in range(0, interval + 1)]
                        p1lst = [epsilon if 0 == p1 else (1 - epsilon) for _ in range(0, interval + 1)]
                        p2lst = [epsilon if 0 == p2 else (1 - epsilon) for _ in range(0, interval + 1)]
                        p3lst = [epsilon if 0 == p3 else (1 - epsilon) for _ in range(0, interval + 1)]
                        if 0 == d:
                            p2list.append([np.array(steplst), np.array(p1lst), np.array(p2lst), np.array(p3lst)])
                        elif 1 == d:
                            p2list.append([np.array(p1lst), np.array(steplst), np.array(p2lst), np.array(p3lst)])
                        elif 2 == d:
                            p2list.append([np.array(p1lst), np.array(p2lst), np.array(steplst), np.array(p3lst)])
                        else:
                            p2list.append([np.array(p1lst), np.array(p2lst), np.array(p3lst), np.array(steplst)])
                    p1list.append(p2list)
                dlist.append(p1list)
            emptyList.append(dlist)
        return emptyList

    @staticmethod
    def __FillLine(func: IntegrandV2,
                   edge: list,
                   xf: complex,
                   xt: complex,
                   yf: complex,
                   yt: complex,
                   zf: complex,
                   zt: complex,
                   wf: complex,
                   wt: complex,
                   interval: int) -> float:
        [x, y, z, w] = np.copy(edge)
        x = xf + (xt - xf) * x
        y = yf + (yt - yf) * y
        z = zf + (zt - zf) * z
        w = wf + (wt - wf) * w
        v = func.vDenom(x, y, z, w)
        u = np.copy(v)
        v = np.delete(v, 0)
        u = np.delete(u, interval)
        return float(np.sum(np.angle(v / u)))

    @staticmethod
    def CheckPolesFast(func: IntegrandV2,
                       edges: list,
                       fromX: complex,
                       toX: complex,
                       fromY: complex,
                       toY: complex,
                       fromZ: complex,
                       toZ: complex,
                       fromW: complex,
                       toW: complex,
                       interval: int = 50) -> bool:
        """
        :return:
        Whether has a pole

        To avoid treat zeros as the pole, please specify the denominator.
        """
        # rescale X
        listArgs = [[[[math.nan for _ in range(0, 2)] for _ in range(0, 2)] for _ in range(0, 2)] for _ in range(0, 4)]
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
                            argLine1 = Integrators4DV2.__FillLine(func,
                                                                  edges[d1][lineArgs1[0]][lineArgs1[1]][lineArgs1[2]],
                                                                  fromX, toX,
                                                                  fromY, toY,
                                                                  fromZ, toZ,
                                                                  fromW, toW,
                                                                  interval)
                            # if maxPhaseChange > maxAllowedPhaseChange:
                            #     return True
                            listArgs[d1][lineArgs1[0]][lineArgs1[1]][lineArgs1[2]] = argLine1
                        else:
                            argLine1 = listArgs[d1][lineArgs1[0]][lineArgs1[1]][lineArgs1[2]]
                        if math.isnan(listArgs[d2][lineArgs2[0]][lineArgs2[1]][lineArgs2[2]]):
                            argLine2 = Integrators4DV2.__FillLine(func,
                                                                  edges[d2][lineArgs2[0]][lineArgs2[1]][lineArgs2[2]],
                                                                  fromX, toX,
                                                                  fromY, toY,
                                                                  fromZ, toZ,
                                                                  fromW, toW,
                                                                  interval)
                            # if maxPhaseChange > maxAllowedPhaseChange:
                            #     return True
                            listArgs[d2][lineArgs2[0]][lineArgs2[1]][lineArgs2[2]] = argLine2
                        else:
                            argLine2 = listArgs[d2][lineArgs2[0]][lineArgs2[1]][lineArgs2[2]]
                        if math.isnan(listArgs[d1][lineArgs3[0]][lineArgs3[1]][lineArgs3[2]]):
                            argLine3 = Integrators4DV2.__FillLine(func,
                                                                  edges[d1][lineArgs3[0]][lineArgs3[1]][lineArgs3[2]],
                                                                  fromX, toX,
                                                                  fromY, toY,
                                                                  fromZ, toZ,
                                                                  fromW, toW,
                                                                  interval)
                            # if maxPhaseChange > maxAllowedPhaseChange:
                            #     return True
                            listArgs[d1][lineArgs3[0]][lineArgs3[1]][lineArgs3[2]] = argLine3
                        else:
                            argLine3 = listArgs[d1][lineArgs3[0]][lineArgs3[1]][lineArgs3[2]]
                        if math.isnan(listArgs[d2][lineArgs4[0]][lineArgs4[1]][lineArgs4[2]]):
                            argLine4 = Integrators4DV2.__FillLine(func,
                                                                  edges[d2][lineArgs4[0]][lineArgs4[1]][lineArgs4[2]],
                                                                  fromX, toX,
                                                                  fromY, toY,
                                                                  fromZ, toZ,
                                                                  fromW, toW,
                                                                  interval)
                            # if maxPhaseChange > maxAllowedPhaseChange:
                            #     return True
                            listArgs[d2][lineArgs4[0]][lineArgs4[1]][lineArgs4[2]] = argLine4
                        else:
                            argLine4 = listArgs[d2][lineArgs4[0]][lineArgs4[1]][lineArgs4[2]]
                        faceArg = abs(argLine1 + argLine2 - argLine3 - argLine4)
                        if faceArg > 0.5:
                            return True
        return False


class SparseGridIntegrator4DV2(Integrators4DV2):

    def __init__(self, fileFolder: str = "../_Data/SparseGrid4D/GaussPatterson",
                 maxOrder: int = 12, epsilon=1.0e-4,
                 fastCheckPole: int = 20,
                 logLevel: LogLevel = LogLevel.General):
        self.epsilon = epsilon
        self.logLevel = logLevel
        [pts, startIndex, endIndex, weights] = CachedPointListAndWeightList(fileFolder, maxOrder)
        self.maxOrder = len(weights)
        pointList = []
        weightList = []
        for i in range(0, self.maxOrder):
            ptAtOrderI = []
            weightAtOrderI = []
            # ======= get function values ========
            for pointIndex in range(startIndex[i], endIndex[i]):
                ptAtOrderI.append([pts[pointIndex].x, pts[pointIndex].y, pts[pointIndex].z, pts[pointIndex].w])
            for pointIndex in range(0, endIndex[i]):
                weightAtOrderI.append(weights[i][pointIndex])
            pointList.append(np.array(ptAtOrderI))
            weightList.append(np.array(weightAtOrderI))
        self.points = pointList
        self.weights = weightList
        self.fastCheckPole = fastCheckPole
        if self.fastCheckPole > 0:
            self.edges = Integrators4DV2.GenerateEdgeList(fastCheckPole)

    def Integrate(self, func: IntegrandV2,
                  fromX: complex, toX: complex,
                  fromY: complex, toY: complex,
                  fromZ: complex, toZ: complex,
                  fromW: complex, toW: complex) -> [bool, complex]:
        if self.logLevel >= LogLevel.Verbose:
            print("integrate4d: x: ", fromX, " to ", toX, " y: ", fromY, " to ", toY, " z: ", fromZ, " to ", toZ,
                  " w: ", fromW, " to ", toW)
        if self.fastCheckPole > 0 and func.vDenom is not None:
            if self.CheckPolesFast(func, self.edges, fromX, toX, fromY, toY, fromZ, toZ, fromW, toW, self.fastCheckPole):
                return [False, cmath.nan]
        resOld: complex = 0
        strideX: complex = 0.5 * (toX - fromX)
        strideY: complex = 0.5 * (toY - fromY)
        strideZ: complex = 0.5 * (toZ - fromZ)
        strideW: complex = 0.5 * (toW - fromW)
        delta = 0
        values = np.array([])
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            xV = strideX * (self.points[i][:, 0] + 1) + fromX
            yV = strideY * (self.points[i][:, 1] + 1) + fromY
            zV = strideZ * (self.points[i][:, 2] + 1) + fromZ
            wV = strideW * (self.points[i][:, 3] + 1) + fromW
            newValues = func.vFunc(xV, yV, zV, wV)
            values = np.append(values, newValues)
            resNew = complex(np.dot(values, self.weights[i]))
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
                "Sparse Grid 4d failed due to max iteration reached\n, x:{} to {} y:{} to {} z:{} to {} w:{} to {}  last value is "
                .format(fromX, toX, fromY, toY, fromZ, toZ, fromW, toW),
                strideX * strideY * strideZ * strideW * resOld,
                " delta = {}/{}".format(delta, self.epsilon))
        return [False, strideX * strideY * strideZ * strideW * resOld]

    def Finish(self):
        return
