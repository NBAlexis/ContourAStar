import cmath
from enum import IntEnum

from Contour1D.CGeneralGrid import CGeneralGrid
from Contour1D.CGrids import GridDir, GridNeighbour, GridState, OneGrid, AStarResult
from Contour1D.CommonDefinitions import LogLevel
from Contour2D import Integrator2D
from Contour2D.Integrand2D import Integrand2D
from Contour3D import Integrator3D
from Contour3D.Integrand3D import Integrand3D


class CMarkedGrids3D:
    """
    确定y,z path
    要求1：对于x=1,-1,y=1,-1，z能积分。
    要求2: 对于x=1,-1， y,z能2重积分。
    第一步： 对于x=1,-1,y=1,-1,用1维积分确定z的路径。
    第二步： 对x=1, -1，确定的z路径，用2维积分确定y的路径。
    第三步： 对确定的y,z路径，用3维积分确定x的路径
    如果失败，回到第二步，标记一个disconnect，重新寻路，直至失败。

    回到第一步，编辑一个disconnect
    """

    """
    假设width, height都是大于等于3的整数
    假设edge < (width - 1) / 2
    """

    def __init__(self, width: int, height: int, edge: int,
                 integrator: Integrator3D, integrand: Integrand3D, requireEdge: bool = True,
                 maxStep: int = 100000, logLevel: LogLevel = LogLevel.General):
        self.integrator = integrator
        self.integrand = integrand
        self.XGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.YGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.ZGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        vLeft = self.integrator.GetLeftEdgeX()
        self.edgePoints = [vLeft, -vLeft]
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.AStarRes = AStarResult.Unknown
        self.resV = 0
        self.requireEdge = requireEdge

    def Integrate(self) -> [AStarResult, complex]:
        """
        第一步： 对于x=1,-1,y=1,-1,用1维积分确定z的路径。
        """
        while True:
            [aStarRes, self.ZPath, _] = self.ZGrid.FindPath(self.CalculateConnectionStep1)
            print("=================================")
            print(aStarRes, self.ZPath)
            if AStarResult.Finished == aStarRes:
                self.YGrid.ResetGrid()
                [step2Res, result] = self.Step2()
                if AStarResult.Finished == step2Res:
                    self.AStarRes = AStarResult.Finished
                    self.resV = result
                    return [AStarResult.Finished, result]
                else:
                    self.ZGrid.EliminateMiddlePointOfPath()
            else:
                return [AStarResult.Failed, cmath.nan]

    def Step2(self) -> [AStarResult, complex]:
        """
        第二步： 对x=1, -1，确定的z路径，用2维积分确定y的路径。
        """
        while True:
            [aStarRes, self.YPath, _] = self.YGrid.FindPath(self.CalculateConnectionStep2)
            if AStarResult.Finished == aStarRes:
                [step3Res, result] = self.Step3()
                if AStarResult.Finished == step3Res:
                    return [AStarResult.Finished, result]
                else:
                    self.YGrid.EliminateMiddlePointOfPath()
            else:
                return [AStarResult.Failed, cmath.nan]

    def Step3(self) -> [AStarResult, complex]:
        self.XGrid.ResetGrid()
        [aStarRes, self.XPath, res] = self.XGrid.FindPath(self.CalculateConnectionStep3)
        return [aStarRes, res]

    def SetIntegrand(self, integrand: Integrand3D):
        self.integrand = integrand

    def CalculateConnectionStep1(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        第一步： 对于x=1,-1,y=1,-1,用1维积分确定z的路径。
        """
        if not self.requireEdge:
            return [True, 0]
        bHasValue: bool = True
        fromV: complex = gridFrom.v
        toV: complex = gridTo.v
        for xValue in self.edgePoints:
            for yValue in self.edgePoints:
                [bHasValue1, _] = self.integrator.IntegrateZ(
                    self.integrand, xValue, yValue, fromV, toV)
                bHasValue = bHasValue and bHasValue1
                if not bHasValue:
                    break
            if not bHasValue:
                break
        return [bHasValue, 0]

    def CalculateConnectionStep2(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        第二步： 对x=1, -1，确定的z路径，用2维积分确定y的路径。
        """
        if not self.requireEdge:
            return [True, 0]
        bHasValue: bool = True
        fromV: complex = gridFrom.v
        toV: complex = gridTo.v
        for xValue in self.edgePoints:
            for i in range(0, len(self.ZPath) - 1):
                [bHasValue1, _] = self.integrator.IntegrateYZ(self.integrand,
                                                              xValue,
                                                              fromV,
                                                              toV,
                                                              self.ZPath[i],
                                                              self.ZPath[i + 1])
                bHasValue = bHasValue and bHasValue1
                if not bHasValue:
                    break
            if not bHasValue:
                break
        return [bHasValue, 0]

    def CalculateConnectionStep3(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        bHasValue: bool = True
        value: complex = 0
        fromV: complex = gridFrom.v
        toV: complex = gridTo.v
        for i in range(0, len(self.YPath) - 1):
            for j in range(0, len(self.ZPath) - 1):
                [bHasValue1, valueCube] = self.integrator.Integrate(
                    self.integrand,
                    fromV,
                    toV,
                    self.YPath[i], self.YPath[i + 1],
                    self.ZPath[j], self.ZPath[j + 1])
                bHasValue = bHasValue and bHasValue1
                value = value + valueCube
                if not bHasValue:
                    break
            if not bHasValue:
                break
        return [bHasValue, value]

    def GatherInfo(self) -> str:
        """
        assume the integrate is finished
        """
        if AStarResult.Finished == self.AStarRes:
            res = ""
            listPlotX = "zr={"
            listPlotY = "zi={"
            for i in range(0, len(self.ZPath)):
                listPlotX = listPlotX + (
                    ",{}".format(self.ZPath[i].real) if 0 != i else "{}".format(self.ZPath[i].real))
                listPlotY = listPlotY + (
                    ",{}".format(self.ZPath[i].imag) if 0 != i else "{}".format(self.ZPath[i].imag))
            listPlotX = listPlotX + "};\n"
            listPlotY = listPlotY + "};\n"
            res = res + listPlotX + listPlotY
            listPlotX = "yr={"
            listPlotY = "yi={"
            for i in range(0, len(self.YPath)):
                listPlotX = listPlotX + (
                    ",{}".format(self.YPath[i].real) if 0 != i else "{}".format(self.YPath[i].real))
                listPlotY = listPlotY + (
                    ",{}".format(self.YPath[i].imag) if 0 != i else "{}".format(self.YPath[i].imag))
            listPlotX = listPlotX + "};\n"
            listPlotY = listPlotY + "};\n"
            res = res + listPlotX + listPlotY
            listPlotX = "xr={"
            listPlotY = "xi={"
            for i in range(0, len(self.XPath)):
                listPlotX = listPlotX + (
                    ",{}".format(self.XPath[i].real) if 0 != i else "{}".format(self.XPath[i].real))
                listPlotY = listPlotY + (
                    ",{}".format(self.XPath[i].imag) if 0 != i else "{}".format(self.XPath[i].imag))
            listPlotX = listPlotX + "};\n"
            listPlotY = listPlotY + "};\n"
            res = res + listPlotX + listPlotY
            res = res + "Print[\"expecting: {}\"]".format(self.resV)
            res = res + """
res = Sum[NIntegrate[g[x, y, z], 
{x, xr[[u]] + xi[[u]] I, xr[[u + 1]] + xi[[u + 1]] I}, 
{y, yr[[v]] + yi[[v]] I, yr[[v + 1]] + yi[[v + 1]] I}, 
{z, zr[[w]] + zi[[w]] I, zr[[w + 1]] + zi[[w + 1]] I}], 
{u, 1, Length[xr] - 1}, {v, 1, Length[yr] - 1}, {w, 1, Length[zr] - 1}]
"""
            res = res + "ListLinePlot[Transpose[{xr, xi}], AxesLabel -> {\"Re[x]\", \"Im[x]\"}]\n"
            res = res + "ListLinePlot[Transpose[{yr, yi}], AxesLabel -> {\"Re[y]\", \"Im[y]\"}]\n"
            res = res + "ListLinePlot[Transpose[{zr, zi}], AxesLabel -> {\"Re[z]\", \"Im[z]\"}]\n"
            return "(* =========== Copy these to Mathematica ========== *)\n\n" + self.integrand.GetDebugInfo() + "\n" + res
        return "Integrating: {}\n".format(self.integrand.GetDebugInfo())

