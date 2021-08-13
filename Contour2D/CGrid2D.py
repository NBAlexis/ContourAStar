import cmath
from enum import IntEnum

from Contour1D.CGrids import GridDir, GridNeighbour, GridState, OneGrid, AStarResult
from Contour1D.CommonDefinitions import LogLevel
from Contour2D import Integrator2D
from Contour2D.Integrand2D import Integrand2D


class CGrids2D:
    """
    假设width, height都是大于等于3的整数
    假设edge < (width - 1) / 2
    """

    def __init__(self, width: int, height: int, edge: int, integrator: Integrator2D, integrand: Integrand2D,
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
        self.YPath = []
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

    def SetIntegrand(self, integrand: Integrand2D):
        self.integrand = integrand

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
            self.YPath = points
            return [res, points]
        self.lastError = "A path for integrate g[-1, y], or g[1, y] not found."
        return [AStarResult.Failed, []]

    def OneStep2(self) -> AStarResult:
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
                self.CalculateConnectionStep2(i, smallestFnNode.x, smallestFnNode.y)
            if smallestFnNode.neighbourState[i] == GridNeighbour.Connected:
                neighbour = smallestFnNode.neighbourNode[i]
                if neighbour.state != GridState.ClosedList:
                    neighbour.UpdateFn(smallestFnNode, i)
        smallestFnNode.state = GridState.ClosedList
        return AStarResult.Unknown

    def FindPathStep2(self) -> AStarResult:
        self.ResetGrid()
        res: AStarResult = self.OneStep2()
        while AStarResult.Unknown == res and self.steps < self.maxStep:
            self.steps = self.steps + 1
            if self.logLevel >= LogLevel.Verbose:
                print("Step:{}".format(self.steps))
            res = self.OneStep2()
        return res

    def IntegrateStep2(self) -> [AStarResult, complex]:
        res: AStarResult = self.FindPathStep2()
        resV: complex = 0
        if AStarResult.Finished == res:
            endNote = self.gridArray[self.endY][self.endX]
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                resV = resV + endNote.neighbourV[endNote.parentDir]
                endNote = endNote.parent
            return [res, -resV]
        self.lastError = "A path to integrate g[x, y] not found."
        return [AStarResult.Failed, cmath.nan]

    def Integrate(self) -> [AStarResult, complex]:
        maxTry = len(self.allXPoints)
        for i in range(0, maxTry):
            print("========= try {} ==========".format(i))
            self.ResetGrid()
            self.XPoints = self.GetEdgeXPointList(i)
            self.YPath = []
            [res, _] = self.IntegrateStep1()
            if res != AStarResult.Finished:
                continue
            [res, resV] = self.IntegrateStep2()
            if res == AStarResult.Finished:
                return [res, resV]
        return [AStarResult.Failed, cmath.nan]

    def CalculateConnectionStep1(self, direction: GridDir, x: int, y: int):
        if self.gridArray[y][x].neighbourNode[direction] is None:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            return
        if self.gridArray[y][x].neighbourState[direction] == GridNeighbour.Connected:
            return
        oppositeD: GridDir = GridDir.Down - direction
        bHasValue = True
        for xValue in self.XPoints:
            [bHasValue1, _] = self.integrator.PartialIntegrateY(
                self.integrand, xValue, self.gridArray[y][x].v, self.gridArray[y][x].neighbourNode[direction].v)
            bHasValue = bHasValue and bHasValue1
            if not bHasValue:
                break
        if bHasValue:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourV[direction] = 0
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourNode[direction].neighbourV[oppositeD] = -0
        else:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

    def CalculateConnectionStep2(self, direction: GridDir, x: int, y: int):
        if self.gridArray[y][x].neighbourNode[direction] is None:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            return
        if self.gridArray[y][x].neighbourState[direction] == GridNeighbour.Connected:
            return
        oppositeD: GridDir = GridDir.Down - direction
        bHasValue: bool = True
        value: complex = 0
        for i in range(0, len(self.YPath) - 1):
            [bHasValueSurface, valueSurface] = self.integrator.Integrate(self.integrand, self.gridArray[y][x].v,
                                                                         self.gridArray[y][x].neighbourNode[
                                                                             direction].v, self.YPath[i],
                                                                         self.YPath[i + 1])
            bHasValue = bHasValue and bHasValueSurface
            value = value + valueSurface
            if x == self.startX and y == self.startY and direction == GridDir.Down:
                print(valueSurface)
        if bHasValue:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourV[direction] = value
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourNode[direction].neighbourV[oppositeD] = -value
        else:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

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
            endNote = self.gridArray[self.endY][self.endX]
            resV = 0
            finalV = 0
            lastDir = endNote.parentDir
            points = [endNote.v]
            values = []
            while (endNote is not None) and (not endNote.start):
                if lastDir != endNote.parentDir:
                    points.insert(0, endNote.v)
                    lastDir = endNote.parentDir
                    values.insert(0, resV)
                    finalV = finalV + resV
                    resV = 0
                resV = resV + endNote.neighbourV[endNote.parentDir]
                endNote = endNote.parent
            finalV = finalV + resV
            values.insert(0, resV)
            points.insert(0, endNote.v)
            strRes = "res ="
            intRes = ""
            for i in range(0, len(values)):
                for j in range(0, len(self.YPath) - 1):
                    intRes = intRes + "{}NIntegrate[g[x, y], {}x, {}, {}{}, {}y, {}, {}{}" \
                        .format("Print[\"Expecting:{}\"]\nres{} = ".format(-values[i], i) if 0 == j else " + ",
                                "{", points[i], points[i + 1], "}", "{", self.YPath[j], self.YPath[j + 1],
                                "}]" if j < len(self.YPath) - 2 else "}]\n")
                if 0 == i:
                    strRes = strRes + " res{}".format(i)
                else:
                    strRes = strRes + " + res{}".format(i)
            intRes = self.integrand.GetDebugInfo() + "\n" + intRes + "Print[\"Final result, expecting:{}\"]\n".format(
                -finalV) + strRes
            listPlotX = ""
            listPlotY = ""
            for i in range(0, len(self.YPath)):
                listPlotX = listPlotX + (
                    ",{}".format(self.YPath[i].real) if 0 != i else "{}".format(self.YPath[i].real))
                listPlotY = listPlotY + (
                    ",{}".format(self.YPath[i].imag) if 0 != i else "{}".format(self.YPath[i].imag))
            listPlotA = "\nListLinePlot[Transpose[{}{}{}{},{}{}{}{}], AxesLabel -> {}\"Re[y]\", \"Im[y]\"{}]\n" \
                .format("{", "{", listPlotX, "}", "{", listPlotY, "}", "}", "{", "}")
            listPlotX = ""
            listPlotY = ""
            for i in range(0, len(points)):
                listPlotX = listPlotX + (",{}".format(points[i].real) if 0 != i else "{}".format(points[i].real))
                listPlotY = listPlotY + (",{}".format(points[i].imag) if 0 != i else "{}".format(points[i].imag))
            listPlotB = "ListLinePlot[Transpose[{}{}{}{},{}{}{}{}], AxesLabel -> {}\"Re[x]\", \"Im[x]\"{}]\n" \
                .format("{", "{", listPlotX, "}", "{", listPlotY, "}", "}", "{", "}")
            intRes = intRes.replace("j", " I")
            return "(* =========== Copy these to Mathematica ========== *)\n\n" + intRes + listPlotA + listPlotB
        return "Integrating: {}\nLast Error: {}".format(self.integrand.GetDebugInfo(), self.lastError)

    """
    Assume the integrate is finished
    """

    def ShowIntegralPath(self):
        if self.AStarRes != AStarResult.Finished:
            print("AStar not finished or succeed!")
            return
        import matplotlib.pyplot as plt
        endNote = self.gridArray[self.endY][self.endX]
        fv: complex = 0
        while (endNote is not None) and (not endNote.start):
            x1 = endNote.v.real
            y1 = endNote.v.imag
            x2 = endNote.parent.v.real
            y2 = endNote.parent.v.imag
            resV = -endNote.neighbourV[endNote.parentDir]
            fv = fv + resV
            plt.plot([x1, x2], [y1, y2], color='blue')
            plt.text(0.5 * (x1 + x2), 0.5 * (y1 + y2),
                     "({0:.4f} {1} {2:.4f}i)".format(resV.real, '+-'[resV.imag < 0], abs(resV.imag)),
                     ha='center', va='center',
                     bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'none', 'pad': 1})
            endNote = endNote.parent
        plt.title("{}".format(fv))
        plt.show()
