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
from Contour1D.CGeneralGridV2 import CGeneralGridV2, CPathPair, CPathV2, CIntervalListV2
from Contour1D.CGrids import AStarResult, OneGrid
from Contour1D.CommonDefinitions import LogLevel
from Contour1D.IntegrandV2 import IntegrandV2
from Contour3D import Integrator3DV2
# from line_profiler_pycharm import profile


class ExtendPath3DV2:

    def __init__(self, width: int, height: int, edge: int, integrator: Integrator3DV2,
                 integrand: IntegrandV2 = None,
                 maxStep: int = 100000,
                 logLevel: LogLevel = LogLevel.General):
        self.logLevel = logLevel
        self.XGrid = CGeneralGridV2(width, height, edge, maxStep, logLevel)
        self.YGrid = CGeneralGridV2(width, height, edge, maxStep, logLevel)
        self.ZGrid = CGeneralGridV2(width, height, edge, maxStep, logLevel)
        gridSize = width * height
        self.XGrid.SetHash(gridSize ** 5, gridSize ** 4)
        self.YGrid.SetHash(gridSize ** 3, gridSize ** 2)
        self.ZGrid.SetHash(gridSize, 1)
        self.totalIntegrate = ((width - 1) * height + width * (height - 1)) ** 3
        self.integrator = integrator
        self.integrand = integrand
        self.resV = 0j
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.XIntervalDic = {}
        self.YIntervalDic = {}
        self.YZPathDic = {}
        self.YZPathList = []
        self.YZPathListIdx = 0
        self.ZPathDic = {}
        self.ZPathList = []
        self.ZPathListIdx = 0
        self.ConsideringXInterval = None
        self.ConsideringYInterval = None
        self.ConsideringZPath = None
        self.ConsideringYZPath = [None, None]
        self.IntegrateDic = {}

    def Integrate(self, integrand: IntegrandV2 = None) -> [AStarResult, complex]:
        if integrand is not None:
            self.integrand = integrand
        if self.integrand is None:
            return [AStarResult.Finished, 0j]
        self.resV = 0j
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.XIntervalDic = {}
        self.YIntervalDic = {}
        self.YZPathDic = {}
        self.YZPathList = []
        self.YZPathListIdx = 0
        self.ZPathDic = {}
        self.ZPathList = []
        self.ZPathListIdx = 0
        self.ConsideringXInterval = None
        self.ConsideringYInterval = None
        self.ConsideringZPath = None
        self.ConsideringYZPath = [None, None]
        self.IntegrateDic = {}
        # =============================
        self.FindNewYZPair()
        while self.YZPathListIdx < len(self.YZPathList):
            aStarRes = self.ConsiderOneYZPath()
            if aStarRes == AStarResult.Finished:
                return [AStarResult.Finished, self.resV]
        return [AStarResult.Failed, 0j]

    def AddYZPath(self, YPath: CPathV2, ZPath: CPathV2, XInterval):
        pathPair = CPathPair([YPath, ZPath])
        if pathPair not in self.YZPathDic:
            self.YZPathDic[pathPair] = pathPair
            self.YZPathList.append([YPath, ZPath, XInterval])

    def AddZPath(self, path: CPathV2, YInterval):
        if path not in self.ZPathDic:
            self.ZPathDic[path] = path
            self.ZPathList.append([path, YInterval])

    def ConsiderOneZPath(self) -> AStarResult:
        """
        持续运行直到self.ZPathListIdx = len(self.ZPathList)
        """
        if self.ZPathListIdx >= len(self.ZPathList):
            return AStarResult.Failed
        self.ConsideringZPath = self.ZPathList[self.ZPathListIdx][0]
        yInterval = self.ZPathList[self.ZPathListIdx][1]
        self.ZPathListIdx = self.ZPathListIdx + 1
        if self.logLevel >= LogLevel.Verbose:
            print("Z List: {0} / {1}".format(self.ZPathListIdx, len(self.ZPathList)))
        self.YGrid.ResetGrid()
        [aStarRes, nodes, points, _] = self.YGrid.FindPathWithNodes(self.CalculateConnectionFindY)
        if aStarRes != AStarResult.Finished:
            extendedIntervals = self.YGrid.GetAllDisconnectIntervals()
            for interval in extendedIntervals:
                newIntervalList = CIntervalListV2([interval]) if (yInterval is None) else yInterval.AddOne(interval)
                if newIntervalList is not None and newIntervalList not in self.YIntervalDic:
                    self.YIntervalDic[newIntervalList] = newIntervalList
                    self.ConsideringYInterval = newIntervalList
                    self.FindNewZPath()
        else:
            # 成品y
            self.AddYZPath(CPathV2(nodes, points), self.ConsideringZPath, self.ConsideringXInterval)
        return aStarRes

    def FindNewZPath(self):
        self.ZGrid.ResetGrid()
        [aStarResZ, nodes, points, _] = self.ZGrid.FindPathWithNodes(self.CalculateConnectionFindZ)
        if aStarResZ == AStarResult.Finished:
            # 成品z
            newZPath = CPathV2(nodes, points)
            self.AddZPath(newZPath, self.ConsideringYInterval)

    def ConsiderOneYZPath(self) -> AStarResult:
        """
        对于已有的y-z对，寻找x，如果找不到的话，用extend来添加新的y-z对
        """
        if self.YZPathListIdx >= len(self.YZPathList):
            return AStarResult.Failed
        self.ConsideringYZPath = [self.YZPathList[self.YZPathListIdx][0], self.YZPathList[self.YZPathListIdx][1]]
        xInterval = self.YZPathList[self.YZPathListIdx][2]
        self.YZPathListIdx = self.YZPathListIdx + 1
        if self.logLevel >= LogLevel.General:
            print("YZ List: {0} / {1}".format(self.YZPathListIdx, len(self.YZPathList)))
        self.XGrid.ResetGrid()
        [aStarRes, _, points, v] = self.XGrid.FindPathWithNodes(self.CalculateConnectionFindX)
        if aStarRes != AStarResult.Finished:
            allIntervals = self.XGrid.GetAllDisconnectIntervals()
            for interval in allIntervals:
                newIntervalList = CIntervalListV2([interval]) if (xInterval is None) else xInterval.AddOne(interval)
                if newIntervalList is not None and newIntervalList not in self.XIntervalDic:
                    self.XIntervalDic[newIntervalList] = newIntervalList
                    self.ConsideringXInterval = newIntervalList
                    self.FindNewYZPair()
        else:
            self.resV = v
            self.XPath = points
            self.YPath = self.ConsideringYZPath[0].points
            self.ZPath = self.ConsideringYZPath[1].points
            self.AStarRes = AStarResult.Finished
        return aStarRes

    def FindNewYZPair(self):
        self.ZPathDic = {}
        self.ZPathListIdx = 0
        self.ZPathList = []
        self.YIntervalDic = {}
        # add first ZPath
        self.ConsideringYInterval = None
        self.FindNewZPath()
        while self.ZPathListIdx < len(self.ZPathList):
            # self.ConsiderOneZPath()
            aStarResult = self.ConsiderOneZPath()
            if AStarResult.Finished == aStarResult:
                break

    def CalculateConnectionFindY(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        我们有成品z，半成品x，找一个成品y
        找到的路径加入到yzpathlist，用来寻找成品x
        """
        consideringX: list = [] if self.ConsideringXInterval is None else self.ConsideringXInterval.intervals
        consideringZ: list = [] if self.ConsideringZPath is None else self.ConsideringZPath.nodes
        for i in range(0, len(consideringX)):
            for j in range(0, len(consideringZ) - 1):
                [bHasRes, _] = self.IntegrateCubic(
                    consideringX[i].node1, consideringX[i].node2,
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
        consideringX: list = [] if self.ConsideringXInterval is None else self.ConsideringXInterval.intervals
        consideringY: list = [] if self.ConsideringYInterval is None else self.ConsideringYInterval.intervals
        for i in range(0, len(consideringX)):
            for j in range(0, len(consideringY)):
                [bHasRes, _] = self.IntegrateCubic(
                    consideringX[i].node1, consideringX[i].node2,
                    consideringY[j].node1, consideringY[j].node2,
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

    # @profile
    def IntegrateCubic(self,
                       nodeX1: OneGrid, nodeX2: OneGrid,
                       nodeY1: OneGrid, nodeY2: OneGrid,
                       nodeZ1: OneGrid, nodeZ2: OneGrid) -> [bool, complex]:
        if nodeX1.idx > nodeX2.idx:
            xf = nodeX1.v
            xt = nodeX2.v
            if nodeY1.idx > nodeY2.idx:
                yf = nodeY1.v
                yt = nodeY2.v
                if nodeZ1.idx > nodeZ2.idx:
                    key = nodeX1.idx1 + nodeX2.idx2 + nodeY1.idx1 + nodeY2.idx2 + nodeZ1.idx1 + nodeZ2.idx2
                    finalProd = 1
                    zf = nodeZ1.v
                    zt = nodeZ2.v
                else:
                    key = nodeX1.idx1 + nodeX2.idx2 + nodeY1.idx1 + nodeY2.idx2 + nodeZ1.idx2 + nodeZ2.idx1
                    finalProd = -1
                    zf = nodeZ2.v
                    zt = nodeZ1.v
            else:
                yf = nodeY2.v
                yt = nodeY1.v
                if nodeZ1.idx > nodeZ2.idx:
                    key = nodeX1.idx1 + nodeX2.idx2 + nodeY1.idx2 + nodeY2.idx1 + nodeZ1.idx1 + nodeZ2.idx2
                    finalProd = -1
                    zf = nodeZ1.v
                    zt = nodeZ2.v
                else:
                    key = nodeX1.idx1 + nodeX2.idx2 + nodeY1.idx2 + nodeY2.idx1 + nodeZ1.idx2 + nodeZ2.idx1
                    finalProd = 1
                    zf = nodeZ2.v
                    zt = nodeZ1.v
        else:
            xf = nodeX2.v
            xt = nodeX1.v
            if nodeY1.idx > nodeY2.idx:
                yf = nodeY1.v
                yt = nodeY2.v
                if nodeZ1.idx > nodeZ2.idx:
                    key = nodeX1.idx2 + nodeX2.idx1 + nodeY1.idx1 + nodeY2.idx2 + nodeZ1.idx1 + nodeZ2.idx2
                    finalProd = -1
                    zf = nodeZ1.v
                    zt = nodeZ2.v
                else:
                    key = nodeX1.idx2 + nodeX2.idx1 + nodeY1.idx1 + nodeY2.idx2 + nodeZ1.idx2 + nodeZ2.idx1
                    finalProd = 1
                    zf = nodeZ2.v
                    zt = nodeZ1.v
            else:
                yf = nodeY2.v
                yt = nodeY1.v
                if nodeZ1.idx > nodeZ2.idx:
                    key = nodeX1.idx2 + nodeX2.idx1 + nodeY1.idx2 + nodeY2.idx1 + nodeZ1.idx1 + nodeZ2.idx2
                    finalProd = 1
                    zf = nodeZ1.v
                    zt = nodeZ2.v
                else:
                    key = nodeX1.idx2 + nodeX2.idx1 + nodeY1.idx2 + nodeY2.idx1 + nodeZ1.idx2 + nodeZ2.idx1
                    finalProd = -1
                    zf = nodeZ2.v
                    zt = nodeZ1.v
        resStored = self.IntegrateDic.get(key)
        if resStored is not None:
            return [resStored[0], finalProd * resStored[1]]
        [bHasValue1, v] = self.integrator.Integrate(self.integrand, xf, xt, yf, yt, zf, zt)
        self.IntegrateDic[key] = [bHasValue1, v]
        if self.logLevel >= LogLevel.Verbose:
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
        res = res + "Print[\"expecting: {}\"];\n".format(self.resV)
        res = res + "res = Sum[NIntegrate[" + self.integrand.sIntegrand + """, 
{x, xr[[u]] + xi[[u]] I, xr[[u + 1]] + xi[[u + 1]] I}, 
{y, yr[[v]] + yi[[v]] I, yr[[v + 1]] + yi[[v + 1]] I}, 
{z, zr[[w]] + zi[[w]] I, zr[[w + 1]] + zi[[w + 1]] I}], 
{u, 1, Length[xr] - 1}, {v, 1, Length[yr] - 1}, {w, 1, Length[zr] - 1}]
"""
        res = res + "ListLinePlot[Transpose[{xr, xi}], AxesLabel -> {\"Re[x]\", \"Im[x]\"}, PlotRange -> All]\n"
        res = res + "ListLinePlot[Transpose[{yr, yi}], AxesLabel -> {\"Re[y]\", \"Im[y]\"}, PlotRange -> All]\n"
        res = res + "ListLinePlot[Transpose[{zr, zi}], AxesLabel -> {\"Re[z]\", \"Im[z]\"}, PlotRange -> All]\n"
        return "(* =========== Copy these to Mathematica ========== *)\n\nf[x_,y_,z_]:=" + self.integrand.sFunc + ";\n" + res
