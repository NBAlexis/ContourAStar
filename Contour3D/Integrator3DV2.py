import cmath
import math
import os

import numpy as np

from Contour1D.CommonDefinitions import LogLevel
from Contour1D.IntegrandV2 import IntegrandV2
from SparseGridIntegrators.GaussianPatterson import GaussPatterson
from SparseGridIntegrators.NestedQuadrature import NestedQuadrature
from SparseGridIntegrators.SparseGridCache import CachedPointListAndWeightList, CachedPointListAndWeightList3D
from SparseGridIntegrators.SparseGridGenerator4D import SparseGrid4D


class Integrators3DV2:

    def Integrate(self,
                  func: IntegrandV2,
                  fromX: complex,
                  toX: complex,
                  fromY: complex,
                  toY: complex,
                  fromZ: complex,
                  toZ: complex) -> [bool, complex]:
        pass

    def Finish(self):
        pass

    @staticmethod
    def GenerateEdgeList(interval: int, epsilon: float = 1.0e-9) -> list:
        step = (1 - 2 * epsilon) / interval
        steplst = [epsilon + i * step for i in range(0, interval + 1)]
        zerolst = [epsilon for _ in range(0, interval + 1)]
        onelst = [(1 - epsilon) for _ in range(0, interval + 1)]
        # 1
        yf_zf_X = [np.array(steplst), np.array(zerolst), np.array(zerolst)]
        xt_zf_Y = [np.array(onelst), np.array(steplst), np.array(zerolst)]
        yt_zf_X = [np.array(steplst), np.array(onelst), np.array(zerolst)]
        xf_zf_Y = [np.array(zerolst), np.array(steplst), np.array(zerolst)]
        # 2
        xt_yf_Z = [np.array(onelst), np.array(zerolst), np.array(steplst)]
        yf_zt_X = [np.array(steplst), np.array(zerolst), np.array(onelst)]
        xf_yf_Z = [np.array(zerolst), np.array(zerolst), np.array(steplst)]
        # 3
        xt_yt_Z = [np.array(onelst), np.array(onelst), np.array(steplst)]
        xt_zt_Y = [np.array(onelst), np.array(steplst), np.array(onelst)]
        # 4
        yt_zt_X = [np.array(steplst), np.array(onelst), np.array(onelst)]
        xf_yt_Z = [np.array(zerolst), np.array(onelst), np.array(steplst)]
        # 5
        xf_zt_Y = [np.array(zerolst), np.array(steplst), np.array(onelst)]
        return [yf_zf_X, xt_zf_Y, yt_zf_X, xf_zf_Y,
                xt_yf_Z, yf_zt_X, xf_yf_Z,
                xt_yt_Z, xt_zt_Y,
                yt_zt_X, xf_yt_Z,
                xf_zt_Y]

    @staticmethod
    def __FillLine(func: IntegrandV2,
                   edge: list,
                   xf: complex,
                   xt: complex,
                   yf: complex,
                   yt: complex,
                   zf: complex,
                   zt: complex,
                   interval: int) -> float:
        [x, y, z] = edge
        x = xf + (xt - xf) * x
        y = yf + (yt - yf) * y
        z = zf + (zt - zf) * z
        v = func.vDenom(x, y, z)
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
                       interval: int) -> bool:
        arg_yf_zf_X = Integrators3DV2.__FillLine(func, edges[0], fromX, toX, fromY, toY, fromZ, toZ, interval)
        arg_xt_zf_Y = Integrators3DV2.__FillLine(func, edges[1], fromX, toX, fromY, toY, fromZ, toZ, interval)
        arg_yt_zf_X = Integrators3DV2.__FillLine(func, edges[2], fromX, toX, fromY, toY, fromZ, toZ, interval)
        arg_xf_zf_Y = Integrators3DV2.__FillLine(func, edges[3], fromX, toX, fromY, toY, fromZ, toZ, interval)
        if abs(arg_yf_zf_X + arg_xt_zf_Y - arg_yt_zf_X - arg_xf_zf_Y) > 0.5:
            return True
        arg_xt_yf_Z = Integrators3DV2.__FillLine(func, edges[4], fromX, toX, fromY, toY, fromZ, toZ, interval)
        arg_yf_zt_X = Integrators3DV2.__FillLine(func, edges[5], fromX, toX, fromY, toY, fromZ, toZ, interval)
        arg_xf_yf_Z = Integrators3DV2.__FillLine(func, edges[6], fromX, toX, fromY, toY, fromZ, toZ, interval)
        if abs(arg_yf_zf_X + arg_xt_yf_Z - arg_yf_zt_X - arg_xf_yf_Z) > 0.5:
            return True
        arg_xt_yt_Z = Integrators3DV2.__FillLine(func, edges[7], fromX, toX, fromY, toY, fromZ, toZ, interval)
        arg_xt_zt_Y = Integrators3DV2.__FillLine(func, edges[8], fromX, toX, fromY, toY, fromZ, toZ, interval)
        if abs(arg_xt_zf_Y + arg_xt_yt_Z - arg_xt_zt_Y - arg_xt_yf_Z) > 0.5:
            return True
        arg_yt_zt_X = Integrators3DV2.__FillLine(func, edges[9], fromX, toX, fromY, toY, fromZ, toZ, interval)
        arg_xf_yt_Z = Integrators3DV2.__FillLine(func, edges[10], fromX, toX, fromY, toY, fromZ, toZ, interval)
        if abs(arg_yt_zf_X + arg_xt_yt_Z - arg_yt_zt_X - arg_xf_yt_Z) > 0.5:
            return True
        arg_xf_zt_Y = Integrators3DV2.__FillLine(func, edges[11], fromX, toX, fromY, toY, fromZ, toZ, interval)
        if abs(arg_yf_zt_X + arg_xt_zt_Y - arg_yt_zt_X - arg_xf_zt_Y) > 0.5:
            return True
        if abs(arg_xf_zf_Y + arg_xf_yt_Z - arg_xf_zt_Y - arg_xf_yf_Z) > 0.5:
            return True
        return False


class SparseGridIntegrator3DV2(Integrators3DV2):

    def __init__(self, fileFolder: str = "../_Data/SparseGrid3D/GaussPatterson",
                 maxOrder: int = 12, epsilon=1.0e-4,
                 fastCheckPole: int = 20,
                 logLevel: LogLevel = LogLevel.General):
        self.epsilon = epsilon
        self.logLevel = logLevel
        [pts, startIndex, endIndex, weights] = CachedPointListAndWeightList3D(fileFolder, maxOrder)
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
            self.edges = Integrators3DV2.GenerateEdgeList(fastCheckPole)

    def Integrate(self, func: IntegrandV2,
                  fromX: complex, toX: complex,
                  fromY: complex, toY: complex,
                  fromZ: complex, toZ: complex) -> [bool, complex]:
        if self.logLevel >= LogLevel.Verbose:
            print("integrate3dV2: x: ", fromX, " to ", toX, " y: ", fromY, " to ", toY, " z: ", fromZ, " to ", toZ)
        if self.fastCheckPole > 0 and func.vDenom is not None:
            if self.CheckPolesFast(func, self.edges, fromX, toX, fromY, toY, fromZ, toZ, self.fastCheckPole):
                return [False, cmath.nan]
        resOld: complex = 0
        strideX: complex = 0.5 * (toX - fromX)
        strideY: complex = 0.5 * (toY - fromY)
        strideZ: complex = 0.5 * (toZ - fromZ)
        delta = 0
        values = np.array([])
        for i in range(0, self.maxOrder):
            # ======= get function values ========
            xV = strideX * (self.points[i][:, 0] + 1) + fromX
            yV = strideY * (self.points[i][:, 1] + 1) + fromY
            zV = strideZ * (self.points[i][:, 2] + 1) + fromZ
            newValues = func.vFunc(xV, yV, zV)
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
                        print("Sparse Grid 3d done, value is ", strideX * strideY * strideZ * resNew)
                    return [True, strideX * strideY * strideZ * resNew]
            resOld = resNew
        if self.logLevel >= LogLevel.General:
            print(
                "Sparse Grid 3d failed due to max iteration reached\n, x:{} to {} y:{} to {} z:{} to {} last value is "
                .format(fromX, toX, fromY, toY, fromZ, toZ),
                strideX * strideY * strideZ * resOld,
                " delta = {}/{}".format(delta, self.epsilon))
        return [False, strideX * strideY * strideZ * resOld]

    def Finish(self):
        return
