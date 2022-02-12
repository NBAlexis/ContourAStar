"""
使用:

（1）z=-1到1，y=-1到1 w=-1到1，寻找x的路径。
（2）如果找不到，针对每个extend x，寻找y,z,w的路径。
    （2.1) z,w=-1到1，x intervals，寻找y的路径。
    （2.2）如果没找到，形成extend y。
        （2.2.1） 对w=-1到1， x-y intervals，寻找z的路径
        （2.2.2） 如果没找到，形成extend z。
            （2.2.2.1） x-y-w intervals，寻找能使x-y-w intervals积分的完整的w路径
            （2.2.2.1） 利用找到的w路径，重复(2.2.1)直到找到z的路径
        （2.2.3） 找到能使x-y intervals积分的完整的z和w路径重复(2.1)，直到找到y的路径
    （2.3）找到能使x intervals积分的完整的y,z和w路径，重复重复（1）
（3）如果yzw路径列表为空，则失败。

"""
from Contour1D.CGeneralGridV2 import CPathV2, CPathPair, CGeneralGridV2, CIntervalListV2
from Contour1D.CGrids import AStarResult, OneGrid
from Contour1D.CommonDefinitions import LogLevel
from Contour1D.IntegrandV2 import IntegrandV2
from Contour4D import Integrator4DV2


class ExtendPath4D:

    def __init__(self, width: int, height: int, edge: int, integrator: Integrator4DV2,
                 integrand: IntegrandV2 = None,
                 maxStep: int = 100000,
                 logLevel: LogLevel = LogLevel.General):
        self.logLevel = logLevel
        self.XGrid = CGeneralGridV2(width, height, edge, maxStep, logLevel)
        self.YGrid = CGeneralGridV2(width, height, edge, maxStep, logLevel)
        self.ZGrid = CGeneralGridV2(width, height, edge, maxStep, logLevel)
        self.WGrid = CGeneralGridV2(width, height, edge, maxStep, logLevel)
        gridSize = width * height
        self.totalIntegrate = ((width - 1) * height + width * (height - 1)) ** 4
        self.XGrid.SetHash(gridSize ** 3, gridSize ** 2)
        self.YGrid.SetHash(gridSize, 1)
        self.ZGrid.SetHash(gridSize ** 3, gridSize ** 2)
        self.WGrid.SetHash(gridSize, 1)
        self.integrator = integrator
        self.integrand = integrand
        self.resV = 0j
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.WPath = []
        self.XIntervalDic = {}
        self.YIntervalDic = {}
        self.ZIntervalDic = {}
        self.WPathList = []
        self.ZWPathList = []
        self.YZWPathList = []
        self.WPathDic = {}
        self.ZWPathDic = {}
        self.YZWPathDic = {}
        self.WPathListIdx = 0
        self.ZWPathListIdx = 0
        self.YZWPathListIdx = 0
        self.ConsideringXInterval = None
        self.ConsideringYInterval = None
        self.ConsideringZInterval = None
        self.ConsideringWPath = None
        self.ConsideringZWPath = [None, None]
        self.ConsideringYZWPath = [None, None, None]
        self.IntegrateDic = {}
        self.IntegrateDicCount = 0

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
        self.WPath = []
        self.XIntervalDic = {}
        self.YIntervalDic = {}
        self.ZIntervalDic = {}
        self.WPathList = []
        self.ZWPathList = []
        self.YZWPathList = []
        self.WPathDic = {}
        self.ZWPathDic = {}
        self.YZWPathDic = {}
        self.WPathListIdx = 0
        self.ZWPathListIdx = 0
        self.YZWPathListIdx = 0
        self.ConsideringXInterval = None
        self.ConsideringYInterval = None
        self.ConsideringZInterval = None
        self.ConsideringWPath = None
        self.ConsideringZWPath = [None, None]
        self.ConsideringYZWPath = [None, None, None]
        self.IntegrateDic = {}
        self.IntegrateDicCount = 0
        # =============================
        self.FindNewYZWPair()
        while self.YZWPathListIdx < len(self.YZWPathList):
            aStarRes = self.ConsiderOneYZWPath()
            if aStarRes == AStarResult.Finished:
                return [AStarResult.Finished, self.resV]
        return [AStarResult.Failed, 0j]

    def AddYZWPath(self, YPath: CPathV2, ZPath: CPathV2, WPath: CPathV2, xInterval):
        newPathPair = CPathPair([YPath, ZPath, WPath])
        if newPathPair not in self.YZWPathDic:
            self.YZWPathDic[newPathPair] = newPathPair
            self.YZWPathList.append([YPath, ZPath, WPath, xInterval])

    def AddZWPath(self, ZPath: CPathV2, WPath: CPathV2, yInterval):
        newPathPair = CPathPair([ZPath, WPath])
        if newPathPair not in self.ZWPathDic:
            self.ZWPathDic[newPathPair] = newPathPair
            self.ZWPathList.append([ZPath, WPath, yInterval])

    def AddWPath(self, path: CPathV2, zInterval):
        if path not in self.WPathDic:
            self.WPathDic[path] = path
            self.WPathList.append([path, zInterval])

    def ConsiderOneWPath(self) -> AStarResult:
        """
        进入到这里时，x,y都是半成品，w是成品。
        目标是看看z行不行得通
        """
        if self.WPathListIdx >= len(self.WPathList):
            return AStarResult.Failed
        self.ConsideringWPath = self.WPathList[self.WPathListIdx][0]
        zInterval = self.WPathList[self.WPathListIdx][1]
        self.WPathListIdx = self.WPathListIdx + 1
        if self.logLevel >= LogLevel.General:
            print("        W List: {0} / {1}".format(self.WPathListIdx, len(self.WPathList)))
        self.ZGrid.ResetGrid()
        [aStarRes, nodes, points, _] = self.ZGrid.FindPathWithNodes(self.CalculateConnectionFindZ)
        if aStarRes != AStarResult.Finished:
            intervals = self.ZGrid.GetAllDisconnectIntervals()
            for interval in intervals:
                newIntervalList = CIntervalListV2([interval]) if (zInterval is None) else zInterval.AddOne(interval)
                if newIntervalList is not None and newIntervalList not in self.ZIntervalDic:
                    self.ZIntervalDic[newIntervalList] = newIntervalList
                    self.ConsideringZInterval = newIntervalList
                    self.FindNewWPath()
        else:
            # 成品z
            self.AddZWPath(CPathV2(nodes, points), self.ConsideringWPath, self.ConsideringYInterval)
        return aStarRes

    def FindNewWPath(self):
        self.WGrid.ResetGrid()
        [aStarResZ, nodes, points, _] = self.WGrid.FindPathWithNodes(self.CalculateConnectionFindW)
        if aStarResZ == AStarResult.Finished:
            # 成品w
            newWPath = CPathV2(nodes, points)
            self.AddWPath(newWPath, self.ConsideringZInterval)

    def ConsiderOneYZWPath(self) -> AStarResult:
        """
        对于已有的y-z对，寻找x，如果找不到的话，用extend来添加新的y-z对
        """
        if self.YZWPathListIdx >= len(self.YZWPathList):
            return AStarResult.Failed
        self.ConsideringYZWPath = [self.YZWPathList[self.YZWPathListIdx][0],
                                   self.YZWPathList[self.YZWPathListIdx][1],
                                   self.YZWPathList[self.YZWPathListIdx][2]]
        xInterval = self.YZWPathList[self.YZWPathListIdx][3]
        self.YZWPathListIdx = self.YZWPathListIdx + 1
        if self.logLevel >= LogLevel.General:
            print("YZW List: {0} / {1}".format(self.YZWPathListIdx, len(self.YZWPathList)))
        self.XGrid.ResetGrid()
        [aStarRes, _, points, v] = self.XGrid.FindPathWithNodes(self.CalculateConnectionFindX)
        if aStarRes != AStarResult.Finished:
            intervals = self.XGrid.GetAllDisconnectIntervals()
            for interval in intervals:
                newIntervalList = CIntervalListV2([interval]) if (xInterval is None) else xInterval.AddOne(interval)
                if newIntervalList is not None and newIntervalList not in self.XIntervalDic:
                    self.XIntervalDic[newIntervalList] = newIntervalList
                    self.ConsideringXInterval = newIntervalList
                    self.FindNewYZWPair()
        else:
            self.resV = v
            self.XPath = points
            self.YPath = self.ConsideringYZWPath[0].points
            self.ZPath = self.ConsideringYZWPath[1].points
            self.WPath = self.ConsideringYZWPath[2].points
            self.AStarRes = AStarResult.Finished
        return aStarRes

    def FindNewZWPair(self):
        """
        这个函数的目标：对于已有的 半成品x-path, y-path
        寻找能使 半成品x-path, y-path 积分成功的 z-w path对。
        """
        self.WPathListIdx = 0
        self.WPathList = []
        self.WPathDic = {}
        self.ZIntervalDic = {}
        # add first WPath
        self.ConsideringZInterval = None
        self.FindNewWPath()
        while self.WPathListIdx < len(self.WPathList):
            # self.ConsiderOneWPath()
            aStarResult = self.ConsiderOneWPath()
            if AStarResult.Finished == aStarResult:
                break

    def ConsiderOneZWPath(self) -> AStarResult:
        """
        是对已找到的ZW-path，看看y行不行得通
        """
        if self.ZWPathListIdx >= len(self.ZWPathList):
            return AStarResult.Failed
        self.ConsideringZWPath = [self.ZWPathList[self.ZWPathListIdx][0], self.ZWPathList[self.ZWPathListIdx][1]]
        yInterval = self.ZWPathList[self.ZWPathListIdx][2]
        self.ZWPathListIdx = self.ZWPathListIdx + 1
        if self.logLevel >= LogLevel.General:
            print("    ZW List: {0} / {1}".format(self.ZWPathListIdx, len(self.ZWPathList)))
        self.YGrid.ResetGrid()
        [aStarRes, nodes, points, _] = self.YGrid.FindPathWithNodes(self.CalculateConnectionFindY)
        if aStarRes != AStarResult.Finished:
            intervals = self.YGrid.GetAllDisconnectIntervals()
            for interval in intervals:
                newIntervalList = CIntervalListV2([interval]) if (yInterval is None) else yInterval.AddOne(interval)
                if newIntervalList is not None and newIntervalList not in self.YIntervalDic:
                    self.YIntervalDic[newIntervalList] = newIntervalList
                    self.ConsideringYInterval = newIntervalList
                    self.FindNewZWPair()
        else:
            # 成品y
            self.AddYZWPath(CPathV2(nodes, points), self.ConsideringZWPath[0], self.ConsideringZWPath[1],
                            self.ConsideringXInterval)
        return aStarRes

    def FindNewYZWPair(self):
        """
        这个函数的目标：对于已有的 半成品x-path。
        寻找能使 半成品x-path 积分成功的 y-z-w path对。
        第一次进来时，半成品x-path为None
        """
        self.ZWPathListIdx = 0
        self.ZWPathList = []
        self.ZWPathDic = {}
        self.YIntervalDic = {}
        # add first ZPath
        self.ConsideringYInterval = None
        self.FindNewZWPair()
        while self.ZWPathListIdx < len(self.ZWPathList):
            # Consider ZW Path, 是对已找到的ZW-path，看看y行不行得通
            # self.ConsiderOneZWPath()
            aStarResult = self.ConsiderOneZWPath()
            if AStarResult.Finished == aStarResult:
                break

    def CalculateConnectionFindW(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        consideringX: list = [] if self.ConsideringXInterval is None else self.ConsideringXInterval.intervals
        consideringY: list = [] if self.ConsideringYInterval is None else self.ConsideringYInterval.intervals
        consideringZ: list = [] if self.ConsideringZInterval is None else self.ConsideringZInterval.intervals
        for i in range(0, len(consideringX)):
            for j in range(0, len(consideringY)):
                for k in range(0, len(consideringZ)):
                    [bHasRes, _] = self.IntegrateHyperCubic(
                        consideringX[i].node1, consideringX[i].node2,
                        consideringY[j].node1, consideringY[j].node2,
                        consideringZ[k].node1, consideringZ[k].node2,
                        gridFrom, gridTo)
                    if not bHasRes:
                        return [False, 0j]
        return [True, 0j]

    def CalculateConnectionFindZ(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        consideringX: list = [] if self.ConsideringXInterval is None else self.ConsideringXInterval.intervals
        consideringY: list = [] if self.ConsideringYInterval is None else self.ConsideringYInterval.intervals
        consideringW: list = self.ConsideringWPath.nodes
        for i in range(0, len(consideringX)):
            for j in range(0, len(consideringY)):
                for k in range(0, len(consideringW) - 1):
                    [bHasRes, _] = self.IntegrateHyperCubic(
                        consideringX[i].node1, consideringX[i].node2,
                        consideringY[j].node1, consideringY[j].node2,
                        gridFrom, gridTo,
                        consideringW[k], consideringW[k + 1])
                    if not bHasRes:
                        return [False, 0j]
        return [True, 0j]

    def CalculateConnectionFindY(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        对于已有的y-z, 找合适的x
        """
        consideringX: list = [] if self.ConsideringXInterval is None else self.ConsideringXInterval.intervals
        consideringZ: list = self.ConsideringZWPath[0].nodes
        consideringW: list = self.ConsideringZWPath[1].nodes
        for i in range(0, len(consideringX)):
            for j in range(0, len(consideringZ) - 1):
                for k in range(0, len(consideringW) - 1):
                    [bHasRes, _] = self.IntegrateHyperCubic(
                        consideringX[i].node1, consideringX[i].node2,
                        gridFrom, gridTo,
                        consideringZ[j], consideringZ[j + 1],
                        consideringW[k], consideringW[k + 1])
                    if not bHasRes:
                        return [False, 0j]
        return [True, 0j]

    def CalculateConnectionFindX(self, gridFrom: OneGrid, gridTo: OneGrid) -> [bool, complex]:
        """
        对于已有的y-z, 找合适的x
        """
        vRes: complex = 0j
        consideringY: list = self.ConsideringYZWPath[0].nodes
        consideringZ: list = self.ConsideringYZWPath[1].nodes
        consideringW: list = self.ConsideringYZWPath[2].nodes
        for i in range(0, len(consideringY) - 1):
            for j in range(0, len(consideringZ) - 1):
                for k in range(0, len(consideringW) - 1):
                    [bHasRes, v] = self.IntegrateHyperCubic(
                        gridFrom, gridTo,
                        consideringY[i], consideringY[i + 1],
                        consideringZ[j], consideringZ[j + 1],
                        consideringW[k], consideringW[k + 1])
                    vRes = vRes + v
                    if not bHasRes:
                        return [False, 0j]
        return [True, vRes]

    def IntegrateHyperCubic(self,
                            nodeX1: OneGrid, nodeX2: OneGrid,
                            nodeY1: OneGrid, nodeY2: OneGrid,
                            nodeZ1: OneGrid, nodeZ2: OneGrid,
                            nodeW1: OneGrid, nodeW2: OneGrid) -> [bool, complex]:
        if nodeX1.idx > nodeX2.idx:
            if nodeY1.idx > nodeY2.idx:
                finalProd = 1
                key1 = nodeX1.idx1 + nodeX2.idx2 + nodeY1.idx1 + nodeY2.idx2
            else:
                finalProd = -1
                key1 = nodeX1.idx1 + nodeX2.idx2 + nodeY1.idx2 + nodeY2.idx1
        else:
            if nodeY1.idx > nodeY2.idx:
                finalProd = -1
                key1 = nodeX1.idx2 + nodeX2.idx1 + nodeY1.idx1 + nodeY2.idx2
            else:
                finalProd = 1
                key1 = nodeX1.idx2 + nodeX2.idx1 + nodeY1.idx2 + nodeY2.idx1
        if nodeZ1.idx > nodeZ2.idx:
            if nodeW1.idx > nodeW2.idx:
                key2 = nodeZ1.idx1 + nodeZ2.idx2 + nodeW1.idx1 + nodeW2.idx2
            else:
                finalProd = -1 * finalProd
                key2 = nodeZ1.idx1 + nodeZ2.idx2 + nodeW1.idx2 + nodeW2.idx1
        else:
            if nodeW1.idx > nodeW2.idx:
                finalProd = -1 * finalProd
                key2 = nodeZ1.idx2 + nodeZ2.idx1 + nodeW1.idx1 + nodeW2.idx2
            else:
                key2 = nodeZ1.idx2 + nodeZ2.idx1 + nodeW1.idx2 + nodeW2.idx1
        ywDic = self.IntegrateDic.get(key1)
        if ywDic is None:
            [bHasValue1, v] = self.integrator.Integrate(
                self.integrand,
                nodeX1.v, nodeX2.v,
                nodeY1.v, nodeY2.v,
                nodeZ1.v, nodeZ2.v,
                nodeW1.v, nodeW2.v)
            ywDic = {key2: [bHasValue1, v * finalProd]}
            self.IntegrateDic[key1] = ywDic
            self.IntegrateDicCount = self.IntegrateDicCount + 1
            if self.logLevel >= LogLevel.General:
                print("Dictionary length = ", self.IntegrateDicCount, " / ", self.totalIntegrate)
            return [bHasValue1, v]
        res = ywDic.get(key2)
        if res is None:
            [bHasValue1, v] = self.integrator.Integrate(
                self.integrand,
                nodeX1.v, nodeX2.v,
                nodeY1.v, nodeY2.v,
                nodeZ1.v, nodeZ2.v,
                nodeW1.v, nodeW2.v)
            ywDic[key2] = [bHasValue1, v * finalProd]
            self.IntegrateDicCount = self.IntegrateDicCount + 1
            if self.logLevel >= LogLevel.General:
                print("Dictionary length = ", self.IntegrateDicCount, " / ", self.totalIntegrate)
            return [bHasValue1, v * finalProd]
        return [res[0], res[1] * finalProd]

    def GatherInfo(self) -> str:
        res = ""
        listPlotX = "wr={"
        listPlotY = "wi={"
        for i in range(0, len(self.WPath)):
            listPlotX = listPlotX + (
                ",{}".format(self.WPath[i].real) if 0 != i else "{}".format(self.WPath[i].real))
            listPlotY = listPlotY + (
                ",{}".format(self.WPath[i].imag) if 0 != i else "{}".format(self.WPath[i].imag))
        listPlotX = listPlotX + "};\n"
        listPlotY = listPlotY + "};\n"
        res = res + listPlotX + listPlotY
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
{z, zr[[s]] + zi[[s]] I, zr[[s + 1]] + zi[[s + 1]] I}, 
{w, wr[[t]] + wi[[t]] I, wr[[t + 1]] + wi[[t + 1]] I}], 
{u, 1, Length[xr] - 1}, {v, 1, Length[yr] - 1}, {s, 1, Length[zr] - 1}, {t, 1, Length[wr] - 1}]
"""
        res = res + "ListLinePlot[Transpose[{xr, xi}], AxesLabel -> {\"Re[x]\", \"Im[x]\"}, PlotRange -> All]\n"
        res = res + "ListLinePlot[Transpose[{yr, yi}], AxesLabel -> {\"Re[y]\", \"Im[y]\"}, PlotRange -> All]\n"
        res = res + "ListLinePlot[Transpose[{zr, zi}], AxesLabel -> {\"Re[z]\", \"Im[z]\"}, PlotRange -> All]\n"
        res = res + "ListLinePlot[Transpose[{wr, wi}], AxesLabel -> {\"Re[w]\", \"Im[w]\"}, PlotRange -> All]\n"
        return "(* =========== Copy these to Mathematica ========== *)\n\n" + "f[x_,y_,z_,w_]:=" + self.integrand.sFunc + ";\n" + res
