"""
使用:

（1）y=-1到1，寻找x的路径。
（2）如果找不到，针对每个blocked x，寻找y的路径。
（3）如果（2）找到的y的路径是新的，那么加入到y路径列表。
（4）如果y路径列表为空，则失败。

"""
from Contour1D.CGeneralGrid import CGeneralGrid, CPath, CIntervalList
from Contour1D.CGrids import AStarResult, OneGrid
from Contour1D.CommonDefinitions import LogLevel
from Contour2D import Integrator2D
from Contour2D.Integrand2D import Integrand2D


class ExtendPath2DV2:

    def __init__(self, width: int, height: int, edge: int, integrator: Integrator2D,
                 integrand: Integrand2D = None,
                 maxStep: int = 100000,
                 logLevel: LogLevel = LogLevel.General):
        self.logLevel = logLevel
        self.XGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.YGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.totalIntegrates = ((width - 1) * height + width * (height - 1)) ** 2
        gridSize = width * height
        self.idxHash = [gridSize ** 3, gridSize ** 2, gridSize, 1]
        self.integrator = integrator
        self.integrand = integrand
        self.resV = 0j
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.YPathNodes = []
        self.YPathListIdx = 0
        self.XIntervalList = []
        self.YPathList = []
        self.ConsideringXInterval = None
        self.ConsideringYPath = None
        self.IntegrateDic = {-1: [False, 0j]}

    def Integrate(self, integrand: Integrand2D = None) -> [AStarResult, complex]:
        if integrand is not None:
            self.integrand = integrand
        if self.integrand is None:
            return [AStarResult.Finished, 0j]
        self.resV = 0j
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.YPathNodes = []
        self.YPathListIdx = 0
        self.XIntervalList = []
        self.YPathList = []
        self.ConsideringXInterval = None
        self.ConsideringYPath = None
        self.IntegrateDic = {-1: [False, 0j]}
        self.YGrid.ResetGrid()
        [aStarRes, nodes, points, _] = self.YGrid.FindPathWithNodes(self.CalculateConnectionStep2)
        if aStarRes == AStarResult.Finished:
            newYPath = CPath(nodes, points)
            self.AddYPath(newYPath, None)
        while len(self.YPathList) > self.YPathListIdx:
            [aStarRes, v] = self.ConsiderOnePath()
            if aStarRes == AStarResult.Finished:
                return [aStarRes, v]
        self.AStarRes = AStarResult.Failed
        return [AStarResult.Failed, 0j]

    def AddYPath(self, path: CPath, xInterval):
        exist: bool = False
        for existPath in self.YPathList:
            if existPath[0].EqualTo(path):
                exist = True
                break
        if not exist:
            self.YPathList.append([path, xInterval])

    def ConsiderOnePath(self) -> [AStarResult, complex]:
        if self.YPathListIdx >= len(self.YPathList):
            return AStarResult.Failed
        self.ConsideringYPath = self.YPathList[self.YPathListIdx][0]
        xIntList = self.YPathList[self.YPathListIdx][1]
        self.YPathListIdx = self.YPathListIdx + 1
        if self.logLevel >= LogLevel.General:
            print("Y path list:", self.YPathListIdx, " / ", len(self.YPathList))
        self.XGrid.ResetGrid()
        [aStarRes, _, points, vRes] = self.XGrid.FindPathWithNodes(self.CalculateConnectionStep1)
        if aStarRes != AStarResult.Finished:
            disconnectX = self.XGrid.GetAllDisconnectIntervals()
            for interval in disconnectX:
                newIntervalList = CIntervalList([interval]) if (xIntList is None) else xIntList.AddOne(interval)
                if newIntervalList is not None:
                    bExist = False
                    for intervalList in self.XIntervalList:
                        if newIntervalList.EqualTo(intervalList):
                            bExist = True
                            break
                    if not bExist:
                        self.XIntervalList.append(newIntervalList)
                        self.ConsideringXInterval = newIntervalList
                        self.YGrid.ResetGrid()
                        [aStarResY, nodes, points, _] = self.YGrid.FindPathWithNodes(self.CalculateConnectionStep2)
                        if aStarResY == AStarResult.Finished:
                            newYPath = CPath(nodes, points)
                            self.AddYPath(newYPath, newIntervalList)
        else:
            self.XPath = points
            self.YPath = self.ConsideringYPath.points
            self.resV = vRes
            self.AStarRes = AStarResult.Finished
        return [aStarRes, vRes]

    def CalculateConnectionStep1(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        vRes: complex = 0j
        consideringY: list = [] if self.ConsideringYPath is None else self.ConsideringYPath.nodes
        for i in range(0, len(consideringY) - 1):
            [bHasRes, v] = self.IntegrateCube(gridFrom, gridTo, consideringY[i], consideringY[i + 1])
            if not bHasRes:
                return [False, 0j]
            vRes = vRes + v
        return [True, vRes]

    def CalculateConnectionStep2(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        consideringX: list = [] if self.ConsideringXInterval is None else self.ConsideringXInterval.intervals
        for i in range(0, len(consideringX)):
            [bHasRes, _] = self.IntegrateCube(consideringX[i].node1, consideringX[i].node2, gridFrom, gridTo)
            if not bHasRes:
                return [False, 0j]
        return [True, 0j]

    def IntegrateCube(self, nodeX1: OneGrid, nodeX2: OneGrid, nodeY1: OneGrid, nodeY2: OneGrid) -> [bool, complex]:
        finalProd: complex = 1.0 + 0.0j
        if nodeX1.idx > nodeX2.idx:
            nodeX1, nodeX2 = nodeX2, nodeX1
            finalProd = finalProd * -1.0
        if nodeY1.idx > nodeY2.idx:
            nodeY1, nodeY2 = nodeY2, nodeY1
            finalProd = finalProd * -1.0
        key = nodeX1.idx * self.idxHash[3] + nodeX2.idx * self.idxHash[2] \
            + nodeY1.idx * self.idxHash[1] + nodeY2.idx * self.idxHash[0]
        if key in self.IntegrateDic:
            [bHasValue1, v] = self.IntegrateDic[key]
            return [bHasValue1, v * finalProd]
        [bHasValue1, v] = self.integrator.Integrate(
            self.integrand,
            nodeX1.v, nodeX2.v,
            nodeY1.v, nodeY2.v)
        self.IntegrateDic[key] = [bHasValue1, v]
        if self.logLevel >= LogLevel.General:
            print("Dictionary length = ", len(self.IntegrateDic), " / ", self.totalIntegrates)
        return [bHasValue1, v * finalProd]

    def GatherInfo(self) -> str:
        """
        assume the integrate is finished
        """
        if AStarResult.Finished == self.AStarRes:
            res = ""
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
            res = res + "\nPrint[\"{0} out of {1} integrations considered\"]".format(len(self.IntegrateDic), self.totalIntegrates)
            res = res + """
res = Sum[NIntegrate[g[x, y], 
{x, xr[[u]] + xi[[u]] I, xr[[u + 1]] + xi[[u + 1]] I}, 
{y, yr[[v]] + yi[[v]] I, yr[[v + 1]] + yi[[v + 1]] I}], 
{u, 1, Length[xr] - 1}, {v, 1, Length[yr] - 1}]
"""
            res = res + "ListLinePlot[Transpose[{xr, xi}], AxesLabel -> {\"Re[x]\", \"Im[x]\"}, PlotRange -> All]\n"
            res = res + "ListLinePlot[Transpose[{yr, yi}], AxesLabel -> {\"Re[y]\", \"Im[y]\"}, PlotRange -> All]\n"
            return "(* =========== Copy these to Mathematica ========== *)\n\n" + self.integrand.GetDebugInfo() + "\n" + res
        return "Integrating: {}\n".format(self.integrand.GetDebugInfo())

