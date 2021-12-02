import cmath
from enum import IntEnum

from Contour1D.CGrids import GridDir, GridNeighbour, GridState, OneGrid, AStarResult
from Contour1D.CommonDefinitions import LogLevel
from Contour2D import Integrator2D
from Contour2D.Integrand2D import Integrand2D
from Contour3D import Integrator3D
from Contour3D.Integrand3D import Integrand3D


class CGrids3D:
    """
    第一步： 对每个x,y 用1维积分确定z的路径。
    第二步： 对每个x，确定的z路径，用2维积分确定y的路径。
    第三步： 对确定的y,z路径，用3维积分确定x的路径
    """

    """
    假设width, height都是大于等于3的整数
    假设edge < (width - 1) / 2
    """

    def __init__(self, width: int, height: int, edge: int,
                 integrator: Integrator3D, integrand: Integrand3D,
                 maxStep: int = 100000, logLevel: LogLevel = LogLevel.General):
        self.integrator = integrator
        self.integrand = integrand
        self.midY = height // 2
        self.sep = 2.0 / (width - 1 - 2 * edge)
        self.width = width
        self.height = height
        self.edge = edge
        self.startX = edge
        self.startY = self.midY
        self.endX = width - 1 - edge
        self.endY = self.midY
        self.gridArray = [
            [OneGrid(self, x, y, x + y * width, -1 + (x - edge) * self.sep + (y - self.midY) * self.sep * 1j)
             for x in range(0, width)] for y
            in range(0, height)]
        self.allXPoints = []
        halfSize = (width - 1) // 2
        for i in range(0, halfSize + 1):
            if halfSize - i == edge:
                self.allXPoints.insert(0, -integrator.GetLeftEdgeX())
            else:
                self.allXPoints.append(i * self.sep)
        self.gridList = []
        xoffsets = [0, -1, 1, 0]
        yoffsets = [1, 0, 0, -1]
        for y in range(0, height):
            for x in range(0, width):
                self.gridArray[y][x].SetHn(self.endX, self.endY)
                self.gridList.append(self.gridArray[y][x])
                # fix up the neighbours
                for d in range(0, GridDir.Max):
                    xoffset: int = xoffsets[d] + x
                    yoffset: int = yoffsets[d] + y
                    if 0 <= xoffset < self.width and 0 <= yoffset < self.height:
                        self.gridArray[y][x].SetNeighbour(d, self.gridArray[yoffset][xoffset])
                    else:
                        self.gridArray[y][x].SetNeighbour(d, None)
        self.steps = 0
        self.maxStep = maxStep
        self.logLevel = logLevel
        self.AStarRes = AStarResult.Unknown
        self.XPath = []
        self.YPath = []
        self.ZPath = []
        self.resV = 0
        self.lastError = ""
        self.XPoints = self.GetEdgeXPointList(0)

    def ResetGrid(self):
        for y in range(0, self.height):
            for x in range(0, self.width):
                self.gridArray[y][x].SetHn(self.endX, self.endY)
                self.gridArray[y][x].parent = None
                self.gridArray[y][x].parentDir = GridDir.Max
                self.gridArray[y][x].fn = 0
                self.gridArray[y][x].state = GridState.NotDecide
                # fix up the neighbours
                for d in range(0, GridDir.Max):
                    if self.gridArray[y][x].neighbourNode[d] is None:
                        self.gridArray[y][x].neighbourState[d] = GridNeighbour.NotConnected
                    else:
                        self.gridArray[y][x].neighbourState[d] = GridNeighbour.Unknown
        self.gridArray[self.startY][self.startX].state = GridState.OpenList
        self.gridArray[self.startY][self.startX].start = True
        self.gridArray[self.endY][self.endX].target = True
        self.steps = 0
        self.AStarRes = AStarResult.Unknown
        # for oneGrid in self.gridList:
        #     print("{} {} : {}".format(oneGrid.x, oneGrid.y, oneGrid.state))
        # print("start:", self.startX, self.startY)

    def SetIntegrand(self, integrand: Integrand3D):
        self.integrand = integrand

    """
    第一步： 对每个x,y 用1维积分确定z的路径。
    """

    def CalculateConnectionStep1(self, direction: GridDir, x: int, y: int):
        if self.gridArray[y][x].neighbourNode[direction] is None:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            return
        if self.gridArray[y][x].neighbourState[direction] == GridNeighbour.Connected:
            return
        oppositeD: GridDir = GridDir.Down - direction
        bHasValue = True
        for xValue in self.XPoints:
            for yValue in self.XPoints:
                [bHasValue1, _] = self.integrator.IntegrateZ(
                    self.integrand, xValue, yValue, self.gridArray[y][x].v,
                    self.gridArray[y][x].neighbourNode[direction].v)
                bHasValue = bHasValue and bHasValue1
                if not bHasValue:
                    break
            if not bHasValue:
                break
        if bHasValue:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourV[direction] = 0
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourNode[direction].neighbourV[oppositeD] = 0
        else:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

    """
    第二步： 对每个x,用2维积分确定y的路径。
    """

    def CalculateConnectionStep2(self, direction: GridDir, x: int, y: int):
        if self.gridArray[y][x].neighbourNode[direction] is None:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            return
        if self.gridArray[y][x].neighbourState[direction] == GridNeighbour.Connected:
            return
        oppositeD: GridDir = GridDir.Down - direction
        bHasValue: bool = True
        first2 = [self.XPoints[0], self.XPoints[1]]
        # for xValue in self.XPoints:
        for xValue in first2:
            for i in range(0, len(self.ZPath) - 1):
                [bHasValue1, _] = self.integrator.IntegrateYZ(self.integrand,
                                                              xValue,
                                                              self.gridArray[y][x].v,
                                                              self.gridArray[y][x].neighbourNode[direction].v,
                                                              self.ZPath[i],
                                                              self.ZPath[i + 1])
                bHasValue = bHasValue and bHasValue1
                if not bHasValue:
                    break
            if not bHasValue:
                break
        if bHasValue:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourV[direction] = 0
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourNode[direction].neighbourV[oppositeD] = 0
        else:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

    def CalculateConnectionStep3(self, direction: GridDir, x: int, y: int):
        if self.gridArray[y][x].neighbourNode[direction] is None:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            return
        if self.gridArray[y][x].neighbourState[direction] == GridNeighbour.Connected:
            return
        oppositeD: GridDir = GridDir.Down - direction
        bHasValue: bool = True
        value: complex = 0
        for i in range(0, len(self.YPath) - 1):
            for j in range(0, len(self.ZPath) - 1):
                [bHasValue1, valueCube] = self.integrator.Integrate(
                    self.integrand,
                    self.gridArray[y][x].v,
                    self.gridArray[y][x].neighbourNode[direction].v,
                    self.YPath[i], self.YPath[i + 1],
                    self.ZPath[j], self.ZPath[j + 1])
                bHasValue = bHasValue and bHasValue1
                value = value + valueCube
                if not bHasValue:
                    break
            if not bHasValue:
                break
        if bHasValue:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourV[direction] = value
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourNode[direction].neighbourV[oppositeD] = -value
        else:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

    def OneStep1(self) -> AStarResult:
        smallestFnNode = None
        smallestFn = -1
        for oneGrid in self.gridList:
            if oneGrid.state == GridState.OpenList:
                if oneGrid.fn < smallestFn or smallestFn < 0:
                    smallestFn = oneGrid.fn
                    smallestFnNode = oneGrid
        if smallestFnNode is None:
            self.AStarRes = AStarResult.Failed
            return AStarResult.Failed
        if smallestFnNode.target:
            self.AStarRes = AStarResult.Finished
            return AStarResult.Finished
        for i in range(0, GridDir.Max):
            if smallestFnNode.neighbourState[i] == GridNeighbour.Unknown:
                self.CalculateConnectionStep1(i, smallestFnNode.x, smallestFnNode.y)
            if smallestFnNode.neighbourState[i] == GridNeighbour.Connected:
                neighbour = smallestFnNode.neighbourNode[i]
                if neighbour.state != GridState.ClosedList:
                    neighbour.UpdateFn(smallestFnNode, i)
        smallestFnNode.state = GridState.ClosedList
        return AStarResult.Unknown

    def FindPathStep1(self) -> AStarResult:
        self.ResetGrid()
        res: AStarResult = self.OneStep1()
        while AStarResult.Unknown == res and self.steps < self.maxStep:
            self.steps = self.steps + 1
            if self.logLevel >= LogLevel.Verbose:
                print("Step:{}".format(self.steps))
            res = self.OneStep1()
        return res

    def IntegrateStep1(self) -> [AStarResult, list]:
        """
        第一步，在x=-1, x=1两个点处，对y从-1到1积分。
        找到一条路径使得这两个积分都存在
        """
        res: AStarResult = self.FindPathStep1()
        if AStarResult.Finished == res:
            lastDir = GridDir.Max
            endNote = self.gridArray[self.endY][self.endX]
            points = []
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                if lastDir != endNote.parentDir:
                    points.insert(0, endNote.v)
                    lastDir = endNote.parentDir
                endNote = endNote.parent
            points.insert(0, endNote.v)
            self.ZPath = points
            return [res, points]
        self.lastError = "A path for integrate g[-1, -1, z], or g[1, 1, z] not found."
        return [AStarResult.Failed, []]

    def OneStep2(self) -> AStarResult:
        smallestFnNode = None
        smallestFn = -1
        for oneGrid in self.gridList:
            if oneGrid.state == GridState.OpenList:
                if oneGrid.fn < smallestFn or smallestFn < 0:
                    smallestFn = oneGrid.fn
                    smallestFnNode = oneGrid
        if smallestFnNode is None:
            self.AStarRes = AStarResult.Failed
            return AStarResult.Failed
        if smallestFnNode.target:
            self.AStarRes = AStarResult.Finished
            return AStarResult.Finished
        for i in range(0, GridDir.Max):
            if smallestFnNode.neighbourState[i] == GridNeighbour.Unknown:
                self.CalculateConnectionStep2(i, smallestFnNode.x, smallestFnNode.y)
            if smallestFnNode.neighbourState[i] == GridNeighbour.Connected:
                neighbour = smallestFnNode.neighbourNode[i]
                if neighbour.state != GridState.ClosedList:
                    neighbour.UpdateFn(smallestFnNode, i)
        smallestFnNode.state = GridState.ClosedList
        return AStarResult.Unknown

    def FindPathStep2(self) -> AStarResult:
        self.ResetGrid()
        res: AStarResult = self.OneStep1()
        while AStarResult.Unknown == res and self.steps < self.maxStep:
            self.steps = self.steps + 1
            if self.logLevel >= LogLevel.Verbose:
                print("Step:{}".format(self.steps))
            res = self.OneStep2()
        return res

    def IntegrateStep2(self) -> [AStarResult, list]:
        """
        第一步，在x=-1, x=1两个点处，对y从-1到1积分。
        找到一条路径使得这两个积分都存在
        """
        res: AStarResult = self.FindPathStep1()
        if AStarResult.Finished == res:
            lastDir = GridDir.Max
            endNote = self.gridArray[self.endY][self.endX]
            points = []
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                if lastDir != endNote.parentDir:
                    points.insert(0, endNote.v)
                    lastDir = endNote.parentDir
                endNote = endNote.parent
            points.insert(0, endNote.v)
            self.YPath = points
            return [res, points]
        self.lastError = "A path for integrate g[-1, y, z], or g[1, y, z] not found."
        return [AStarResult.Failed, []]

    def OneStep3(self) -> AStarResult:
        smallestFnNode = None
        smallestFn = -1
        # print("=========")
        for oneGrid in self.gridList:
            if oneGrid.state == GridState.OpenList:
                # print("open list:", oneGrid.x, oneGrid.y)
                if oneGrid.fn < smallestFn or smallestFn < 0:
                    smallestFn = oneGrid.fn
                    smallestFnNode = oneGrid
        # print("=========")
        if smallestFnNode is None:
            self.AStarRes = AStarResult.Failed
            return AStarResult.Failed
        if smallestFnNode.target:
            self.AStarRes = AStarResult.Finished
            return AStarResult.Finished
        for i in range(0, GridDir.Max):
            if smallestFnNode.neighbourState[i] == GridNeighbour.Unknown:
                self.CalculateConnectionStep3(i, smallestFnNode.x, smallestFnNode.y)
            if smallestFnNode.neighbourState[i] == GridNeighbour.Connected:
                neighbour = smallestFnNode.neighbourNode[i]
                if neighbour.state != GridState.ClosedList:
                    neighbour.UpdateFn(smallestFnNode, i)
        smallestFnNode.state = GridState.ClosedList
        return AStarResult.Unknown

    def FindPathStep3(self) -> AStarResult:
        self.ResetGrid()
        res: AStarResult = self.OneStep3()
        while AStarResult.Unknown == res and self.steps < self.maxStep:
            self.steps = self.steps + 1
            if self.logLevel >= LogLevel.Verbose:
                print("Step:{}".format(self.steps))
            res = self.OneStep3()
        return res

    def IntegrateStep3(self) -> [AStarResult, complex]:
        res: AStarResult = self.FindPathStep3()
        print("================ res ============== ", res)
        resV: complex = 0
        if AStarResult.Finished == res:
            endNote = self.gridArray[self.endY][self.endX]
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                resV = resV + endNote.neighbourV[endNote.parentDir]
                endNote = endNote.parent
            # gather X path
            lastDir = GridDir.Max
            endNote = self.gridArray[self.endY][self.endX]
            points = []
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                if lastDir != endNote.parentDir:
                    points.insert(0, endNote.v)
                    lastDir = endNote.parentDir
                endNote = endNote.parent
            points.insert(0, endNote.v)
            self.XPath = points
            self.resV = -resV
            return [res, -resV]
        self.lastError = "A path to integrate g[x, y, z] not found."
        return [AStarResult.Failed, cmath.nan]

    def Integrate(self) -> [AStarResult, complex]:
        maxTry = len(self.allXPoints)
        for i in range(0, maxTry):
            print("========= try {} ==========".format(i))
            self.ResetGrid()
            self.XPoints = self.GetEdgeXPointList(i)
            self.YPath = []
            self.ZPath = []
            [res, _] = self.IntegrateStep1()
            if res != AStarResult.Finished:
                continue
            [res, _] = self.IntegrateStep2()
            if res != AStarResult.Finished:
                continue
            [res, resV] = self.IntegrateStep3()
            if res == AStarResult.Finished:
                return [res, resV]
        return [AStarResult.Failed, cmath.nan]

    def GetEdgeXPointList(self, tryTime: int) -> list:
        vLeft = self.integrator.GetLeftEdgeX()
        if 0 == tryTime:
            return [-vLeft, vLeft]
        elif 1 == tryTime:
            return [-vLeft, 0, vLeft]
        ret = [-vLeft, 0, vLeft]
        for i in range(2, len(self.allXPoints)):
            ret.append(self.allXPoints[i])
            ret.append(-self.allXPoints[i])
        return ret

    def Show(self):
        import matplotlib.pyplot as plt
        sep = self.sep * 0.25
        xoffsets = [0, -sep, sep, 0]
        yoffsets = [sep, 0, 0, -sep]
        for y in range(0, self.height):
            for x in range(0, self.width):
                v = self.gridArray[y][x].v
                xp = v.real
                yp = v.imag
                colorlist = ['black', 'blue', 'red']
                plt.text(xp, yp, "({},{})".format(x, y), color=colorlist[self.gridArray[y][x].state])
                for d in range(0, GridDir.Max):
                    if self.gridArray[y][x].neighbourState[d] == GridNeighbour.Unknown:
                        plt.arrow(x=xp, y=yp, dx=xoffsets[d], dy=yoffsets[d], width=.02, facecolor='yellow')
                    elif self.gridArray[y][x].neighbourState[d] == GridNeighbour.Connected:
                        if self.gridArray[y][x].parentDir == d:
                            plt.arrow(x=xp, y=yp, dx=xoffsets[d], dy=yoffsets[d], width=.02, facecolor='green')
                        else:
                            plt.arrow(x=xp, y=yp, dx=xoffsets[d], dy=yoffsets[d], width=.02, facecolor='blue')
        plt.show()

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
        return "Integrating: {}\nLast Error: {}".format(self.integrand.GetDebugInfo(), self.lastError)


