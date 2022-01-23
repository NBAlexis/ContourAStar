"""
使用:

（1）z=-1到1，y=-1到1，寻找x的路径。
（2）如果找不到，针对每个extend x，寻找y,z的路径。
    （2.1) z=-1到1，x interval，寻找y的路径。
    （2.2）如果没找到，形成extend y。
    （2.3）对每个x,y的interval，找z的路径。
    （2.4）便利z路径列表，找y的路径，对于没找到的情况，更新blocked y interval。如果没找到，回到（2.2）
（3）如果（2）找到的z的路径是新的，那么加入到z路径列表。
（4）如果z路径列表为空，则失败。

"""
from Contour1D.CGeneralGrid import CGeneralGrid, CPath
from Contour1D.CGrids import AStarResult, OneGrid
from Contour1D.CommonDefinitions import LogLevel
from Contour3D import Integrator3D
from Contour3D.Integrand3D import Integrand3D


class ExtendPath3D:

    def __init__(self, width: int, height: int, edge: int, integrator: Integrator3D,
                 integrand: Integrand3D = None,
                 maxStep: int = 100000,
                 logLevel: LogLevel = LogLevel.General):
        self.logLevel = logLevel
        self.XGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.YGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        self.ZGrid = CGeneralGrid(width, height, edge, maxStep, logLevel)
        gridSize = width * height
        self.totalIntegrate = ((width - 1) * height + width * (height - 1)) ** 3
        self.idxHash = [gridSize ** 5, gridSize ** 4, gridSize ** 3, gridSize ** 2, gridSize, 1]
        self.integrator = integrator
        self.integrand = integrand
        self.resV = 0j
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.YPathNodes = []
        self.XPathList = []
        self.YPathList = []
        self.YZPathList = []
        self.YZPathListIdx = 0
        self.ZPathList = []
        self.ZPathListIdx = 0
        self.ConsideringXPath = None
        self.ConsideringYPath = None
        self.ConsideringZPath = None
        self.ConsideringYZPath = [None, None]
        self.IntegrateDic = {-1: [False, 0j]}

    def Integrate(self, integrand: Integrand3D = None) -> [AStarResult, complex]:
        if integrand is not None:
            self.integrand = integrand
        if self.integrand is None:
            return [AStarResult.Finished, 0j]
        self.resV = 0j
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.YPathNodes = []
        self.XPathList = []
        self.YPathList = []
        self.YZPathList = []
        self.YZPathListIdx = 0
        self.ZPathList = []
        self.ZPathListIdx = 0
        self.ConsideringXPath = None
        self.ConsideringYPath = None
        self.ConsideringZPath = None
        self.ConsideringYZPath = [None, None]
        self.IntegrateDic = {-1: [False, 0j]}
        # =============================
        self.FindNewYZPair()
        while self.YZPathListIdx < len(self.YZPathList):
            aStarRes = self.ConsiderOneYZPath()
            if aStarRes == AStarResult.Finished:
                return [AStarResult.Finished, self.resV]
        return [AStarResult.Failed, 0j]

    def AddYZPath(self, YPath: CPath, ZPath: CPath):
        exist: bool = False
        for existPath in self.YZPathList:
            if existPath[0].EqualTo(YPath) and existPath[1].EqualTo(ZPath):
                exist = True
                break
        if not exist:
            self.YZPathList.append([YPath, ZPath])

    def AddZPath(self, path: CPath):
        exist: bool = False
        for existPath in self.ZPathList:
            if existPath.EqualTo(path):
                exist = True
                break
        if not exist:
            self.ZPathList.append(path)

    def ConsiderOneZPath(self) -> AStarResult:
        """
        持续运行直到self.ZPathListIdx = len(self.ZPathList)
        """
        if self.ZPathListIdx >= len(self.ZPathList):
            return AStarResult.Failed
        self.ConsideringZPath = self.ZPathList[self.ZPathListIdx]
        self.ZPathListIdx = self.ZPathListIdx + 1
        if self.logLevel >= LogLevel.Verbose:
            print("Z List: {0} / {1}".format(self.ZPathListIdx, len(self.ZPathList)))
        self.YGrid.ResetGrid()
        [aStarRes, nodes, points, _] = self.YGrid.FindPathWithNodes(self.CalculateConnectionFindY)
        if aStarRes != AStarResult.Finished:
            extendedPaths = self.YGrid.ExtendAllPath()
            for extended in extendedPaths:
                bExist = False
                for path in self.YPathList:
                    if path.EqualTo(extended):
                        bExist = True
                        break
                if not bExist:
                    self.YPathList.append(extended)
                    self.ConsideringYPath = extended
                    self.FindNewZPath()
        else:
            # 成品y
            self.AddYZPath(CPath(nodes, points), self.ConsideringZPath)
        return aStarRes

    def FindNewZPath(self):
        self.ZGrid.ResetGrid()
        [aStarResZ, nodes, points, _] = self.ZGrid.FindPathWithNodes(self.CalculateConnectionFindZ)
        if aStarResZ == AStarResult.Finished:
            # 成品z
            newZPath = CPath(nodes, points)
            self.AddZPath(newZPath)

    def ConsiderOneYZPath(self) -> AStarResult:
        """
        对于已有的y-z对，寻找x，如果找不到的话，用extend来添加新的y-z对
        """
        if self.YZPathListIdx >= len(self.YZPathList):
            return AStarResult.Failed
        self.ConsideringYZPath = self.YZPathList[self.YZPathListIdx]
        self.YZPathListIdx = self.YZPathListIdx + 1
        if self.logLevel >= LogLevel.General:
            print("YZ List: {0} / {1}".format(self.YZPathListIdx, len(self.YZPathList)))
        self.XGrid.ResetGrid()
        [aStarRes, _, points, v] = self.XGrid.FindPathWithNodes(self.CalculateConnectionFindX)
        if aStarRes != AStarResult.Finished:
            extendedPaths = self.XGrid.ExtendAllPath()
            if self.logLevel >= LogLevel.General:
                print("extended X: ", len(extendedPaths))
            for extended in extendedPaths:
                bExist = False
                for path in self.XPathList:
                    if path.EqualTo(extended):
                        bExist = True
                        break
                if not bExist:
                    self.XPathList.append(extended)
                    self.ConsideringXPath = extended
                    self.FindNewYZPair()
        else:
            self.resV = v
            self.XPath = points
            self.YPath = self.ConsideringYZPath[0].points
            self.ZPath = self.ConsideringYZPath[1].points
            self.AStarRes = AStarResult.Finished
        return aStarRes

    def FindNewYZPair(self):
        self.ZPathListIdx = 0
        self.ZPathList = []
        self.YPathList = []
        # add first ZPath
        self.ConsideringYPath = None
        self.FindNewZPath()
        while self.ZPathListIdx < len(self.ZPathList):
            self.ConsiderOneZPath()
            # aStarResult = self.ConsiderOneZPath()
            # if AStarResult.Finished == aStarResult:
            #     break

    def CalculateConnectionFindY(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        我们有成品z，半成品x，找一个成品y
        找到的路径加入到yzpathlist，用来寻找成品x
        """
        consideringX: list = [] if self.ConsideringXPath is None else self.ConsideringXPath.nodes
        consideringZ: list = [] if self.ConsideringZPath is None else self.ConsideringZPath.nodes
        for i in range(0, len(consideringX) - 1):
            for j in range(0, len(consideringZ) - 1):
                [bHasRes, _] = self.IntegrateCubic(
                    consideringX[i], consideringX[i + 1],
                    gridFrom, gridTo,
                    consideringZ[j], consideringZ[j + 1])
                if not bHasRes:
                    return [False, 0j]
        return [True, 0j]

    def CalculateConnectionFindZ(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        我们由半成品路径x, y，找一个成品z，能对x,y积分
        找到的，加入zpathlist，用来寻找成品的y
        """
        consideringX: list = [] if self.ConsideringXPath is None else self.ConsideringXPath.nodes
        consideringY: list = [] if self.ConsideringYPath is None else self.ConsideringYPath.nodes
        for i in range(0, len(consideringX) - 1):
            for j in range(0, len(consideringY) - 1):
                [bHasRes, _] = self.IntegrateCubic(
                    consideringX[i], consideringX[i + 1],
                    consideringY[j], consideringY[j + 1],
                    gridFrom, gridTo)
                if not bHasRes:
                    return [False, 0j]
        return [True, 0j]

    def CalculateConnectionFindX(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        对于已有的y-z, 找合适的x
        """
        vRes: complex = 0j
        consideringY: list = self.ConsideringYZPath[0].nodes
        consideringZ: list = self.ConsideringYZPath[1].nodes
        for i in range(0, len(consideringY) - 1):
            for j in range(0, len(consideringZ) - 1):
                [bHasRes, v] = self.IntegrateCubic(
                    gridFrom, gridTo,
                    consideringY[i], consideringY[i + 1],
                    consideringZ[j], consideringZ[j + 1])
                vRes = vRes + v
                if not bHasRes:
                    return [False, 0j]
        return [True, vRes]

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
        if self.logLevel >= LogLevel.General:
            print("Dictionary length = ", len(self.IntegrateDic), " / ", self.totalIntegrate)
        return [bHasValue1, v * finalProd]

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
        res = res + "ListLinePlot[Transpose[{xr, xi}], AxesLabel -> {\"Re[x]\", \"Im[x]\"}, PlotRange -> All]\n"
        res = res + "ListLinePlot[Transpose[{yr, yi}], AxesLabel -> {\"Re[y]\", \"Im[y]\"}, PlotRange -> All]\n"
        res = res + "ListLinePlot[Transpose[{zr, zi}], AxesLabel -> {\"Re[z]\", \"Im[z]\"}, PlotRange -> All]\n"
        return "(* =========== Copy these to Mathematica ========== *)\n\n" + self.integrand.GetDebugInfo() + "\n" + res

