"""
我们用最繁琐的方式，只求能算出结果，不管效率。

对于一个Y，能connect用y-path并上这个link有解为前提
"""
from Contour1D.CGeneralGrid import CGeneralGrid, TraceBack
from Contour1D.CGrids import AStarResult, OneGrid
from Contour1D.CommonDefinitions import LogLevel
from Contour3D.Integrator3D import Integrators3D
from Contour3D.Integrand3D import Integrand3D


class NestedAStar3D:

    def __init__(self, width: int, height: int, edge: int,
                 integrator: Integrators3D, integrand: Integrand3D = None,
                 maxStep: int = 100000, logLevel: LogLevel = LogLevel.General):
        self.integrator = integrator
        self.integrand = integrand
        self.XGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.YGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.ZGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.AStarRes = AStarResult.Unknown
        self.resV = 0j
        self.CheckXPath = []
        self.CheckYPath = []
        gridSize = width * height
        self.idxHash = [gridSize ** 5, gridSize ** 4, gridSize ** 3, gridSize ** 2, gridSize, 1]
        self.IntegrateDic = {-1: [False, 0j]}
        self.logLevel = logLevel

    def Integrate(self, integrand: Integrand3D = None) -> AStarResult:
        if integrand is not None:
            self.integrand = integrand
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.AStarRes = AStarResult.Unknown
        self.resV = 0j
        self.CheckXPath = []
        self.CheckYPath = []
        self.IntegrateDic = {-1: [False, 0j]}
        self.XGrid.ResetGrid()
        [aStarRes, self.XPath, _] = self.XGrid.FindPath(self.IsXConnected)
        print(self.XPath)
        print(self.YPath)
        print(self.ZPath)
        if self.logLevel >= LogLevel.Verbose:
            self.XGrid.Show()
        if aStarRes == AStarResult.Finished:
            self.FinallyIntegrate()
        return aStarRes

    def IsXConnected(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        这个是用于y-grid的path-finding的。
        """
        self.CheckXPath = TraceBack(gridFrom) + [[gridFrom, gridTo]]
        self.YPath = []
        self.YGrid.ResetGrid()
        [aStarRes, self.YPath, _] = self.YGrid.FindPath(self.IsYConnected)
        if aStarRes == AStarResult.Finished:
            return [True, 0j]
        return [False, 0j]

    def IsYConnected(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        这个是用于y-grid的path-finding的。
        """
        self.CheckYPath = TraceBack(gridFrom) + [[gridFrom, gridTo]]
        self.ZGrid.ResetGrid()
        [aStarRes, self.ZPath, _] = self.ZGrid.FindPath(self.IsZConnected)
        if aStarRes == AStarResult.Finished:
            return [True, 0j]
        return [False, 0j]

    def IsZConnected(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        这个是用于z-grid的path-finding的。
        对于z，要求对y-path,x-path都能积分。
        """
        bHasValue: bool = True
        for xIntervals in self.CheckXPath:
            for yIntervals in self.CheckYPath:
                [bHasValue1, _] = self.IntegrateCubic(
                    xIntervals[0],
                    xIntervals[1],
                    yIntervals[0],
                    yIntervals[1],
                    gridFrom, gridTo)
                bHasValue = bHasValue and bHasValue1
                if not bHasValue:
                    break
            if not bHasValue:
                break
        return [bHasValue, 0]

    def IntegrateCubic(self,
                       nodeX1: OneGrid, nodeX2: OneGrid,
                       nodeY1: OneGrid, nodeY2: OneGrid,
                       nodeZ1: OneGrid, nodeZ2: OneGrid) -> [bool, complex]:
        finalProd: complex = 1.0 + 0.0j
        if nodeX1.idx > nodeX2.idx:
            nodeX1, nodeX2 = nodeX2, nodeX1
            finalProd = finalProd * -1.0
        if nodeY1.idx > nodeY2.idx:
            nodeY1, nodeY2 = nodeY2, nodeY1
            finalProd = finalProd * -1.0
        if nodeZ1.idx > nodeZ2.idx:
            nodeZ1, nodeZ2 = nodeZ2, nodeZ1
            finalProd = finalProd * -1.0
        key = nodeX1.idx * self.idxHash[0] + nodeX2.idx * self.idxHash[1] \
              + nodeY1.idx * self.idxHash[2] + nodeY2.idx * self.idxHash[3] \
              + nodeZ1.idx * self.idxHash[4] + nodeZ2.idx * self.idxHash[5]
        if key in self.IntegrateDic:
            [bHasValue1, v] = self.IntegrateDic[key]
            return [bHasValue1, v * finalProd]
        [bHasValue1, v] = self.integrator.Integrate(
            self.integrand,
            nodeX1.v, nodeX2.v,
            nodeY1.v, nodeY2.v,
            nodeZ1.v, nodeZ2.v)
        self.IntegrateDic[key] = [bHasValue1, v]
        print("Dictionary length = ", len(self.IntegrateDic))
        return [bHasValue1, v * finalProd]

    def FinallyIntegrate(self):
        xPathAll = self.XGrid.GetPathNodes()
        yPathAll = self.YGrid.GetPathNodes()
        zPathAll = self.ZGrid.GetPathNodes()
        self.resV = 0j
        for i in range(0, len(xPathAll) - 1):
            for j in range(0, len(yPathAll) - 1):
                for k in range(0, len(zPathAll) - 1):
                    oneCubic = self.IntegrateCubic(
                        xPathAll[i], xPathAll[i + 1],
                        yPathAll[j], yPathAll[j + 1],
                        zPathAll[k], zPathAll[k + 1])
                    self.resV = self.resV + oneCubic[1]

    def GatherInfo(self) -> str:
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
