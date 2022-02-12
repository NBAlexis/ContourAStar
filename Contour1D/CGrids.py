import cmath
from enum import IntEnum

import numpy as np

from Contour1D.CommonDefinitions import LogLevel
from Contour1D.Integrand import Integrand
from Contour1D.Integrators import Integrators


class GridState(IntEnum):
    NotDecide = 0
    OpenList = 1
    ClosedList = 2


class GridDir(IntEnum):
    Up = 0
    Left = 1
    Right = 2
    Down = 3
    Max = 4


class GridNeighbour(IntEnum):
    NotConnected = 0
    Unknown = 1
    Connected = 2


class AStarResult(IntEnum):
    Unknown = 0
    Finished = 1
    Failed = 2


class OneGrid:

    def __init__(self, owner, x: int, y: int, idx: int, v: complex):
        self.owner = owner
        self.parent = None
        self.parentDir = GridDir.Max
        self.x = x
        self.y = y
        self.v = v
        self.idx = idx
        self.neighbourV = [0j, 0j, 0j, 0j]
        self.neighbourNode = [None, None, None, None]
        self.neighbourState = [GridNeighbour.Unknown, GridNeighbour.Unknown, GridNeighbour.Unknown,
                               GridNeighbour.Unknown]
        self.state = GridState.NotDecide
        self.fn = 0.0
        self.hn = 0.0
        self.gn = 0.0
        self.start = False
        self.target = False
        self.idx1 = idx
        self.idx2 = idx

    def SetNeighbour(self, direction: GridDir, neighbour):
        if neighbour is None:
            self.neighbourState[direction] = GridNeighbour.NotConnected
            return
        self.neighbourNode[direction] = neighbour

    def CalculateConnection(self, direction: GridDir, integrator: Integrators, integrand: Integrand):
        if self.neighbourNode[direction] is None:
            self.neighbourState[direction] = GridNeighbour.NotConnected
            return
        if self.neighbourState[direction] == GridNeighbour.Connected:
            return
        oppositeD: GridDir = GridDir.Down - direction
        [bHasValue, value] = integrator.Integrate(integrand, self.v, self.neighbourNode[direction].v)
        if bHasValue:
            self.neighbourState[direction] = GridNeighbour.Connected
            self.neighbourV[direction] = value
            self.neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            self.neighbourNode[direction].neighbourV[oppositeD] = -value
        else:
            self.neighbourState[direction] = GridNeighbour.NotConnected
            self.neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

    def SetHn(self, targetX, targetY):
        self.hn = abs(self.x - targetX) + abs(self.y - targetY)

    def UpdateFn(self, parentNode, direction: GridDir):
        oppositeD: GridDir = GridDir.Down - direction
        if self.state == GridState.NotDecide:
            self.parent = parentNode
            self.parentDir = oppositeD
            self.gn = 1.0 + parentNode.gn
            self.fn = self.hn + self.gn
            self.state = GridState.OpenList
            return
        newGn = 1.0 + parentNode.gn
        if newGn < self.gn:
            self.gn = newGn
            self.parent = parentNode
            self.parentDir = oppositeD
            self.fn = self.gn + self.hn

    """
    def __float__(self) -> float:
        if self.state == GridState.OpenList:
            return float(self.fn)
        if self.state == GridState.ClosedList:
            return -10.0
        return -1.0

    def __cmp__(self, other):
        me = float(self)
        you = float(other)
        if me > you:
            return 1
        elif me < you:
            return -1
        return 0
        
    def __lt__(self, other): return self.__cmp__(other) < 0
    def __le__(self, other): return self.__cmp__(other) <= 0
    def __eq__(self, other): return self.__cmp__(other) == 0
    def __ne__(self, other): return self.__cmp__(other) != 0
    def __gt__(self, other): return self.__cmp__(other) > 0
    def __ge__(self, other): return self.__cmp__(other) >= 0
    """

    def __gt__(self, other):
        return self.fn > other.fn

    def __le__(self, other):
        return self.fn <= other.fn


class CGrids:
    """
    假设width, height都是大于等于3的整数
    假设edge < (width - 1) / 2
    """

    def __init__(self, width: int, height: int, edge: int, integrator: Integrators, integrand: Integrand,
                 maxStep: int = 1000, logLevel: LogLevel = LogLevel.General):
        self.integrator = integrator
        self.integrand = integrand
        self.midY = height // 2
        self.sep = 1.0 / (width - 1 - 2 * edge)
        self.width = width
        self.height = height
        self.edge = edge
        self.startX = edge
        self.startY = self.midY
        self.endX = width - 1 - edge
        self.endY = self.midY
        self.gridArray = [
            [OneGrid(self, x, y, x + y * width, (x - edge) * self.sep + (y - self.midY) * self.sep * 1j)
             for x in range(0, width)] for y
            in range(0, height)]
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
        self.gridArray[self.startY][self.startX].state = GridState.OpenList
        self.gridArray[self.startY][self.startX].start = True
        self.gridArray[self.endY][self.endX].target = True
        self.steps = 0
        self.maxStep = maxStep
        self.logLevel = logLevel
        self.AStarRes = AStarResult.Unknown

    def OneStep(self) -> AStarResult:
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
                smallestFnNode.CalculateConnection(i, self.integrator, self.integrand)
            if smallestFnNode.neighbourState[i] == GridNeighbour.Connected:
                neighbour = smallestFnNode.neighbourNode[i]
                if neighbour.state != GridState.ClosedList:
                    neighbour.UpdateFn(smallestFnNode, i)
        smallestFnNode.state = GridState.ClosedList
        return AStarResult.Unknown

    def FindPath(self) -> AStarResult:
        res: AStarResult = self.OneStep()
        while AStarResult.Unknown == res and self.steps < self.maxStep:
            self.steps = self.steps + 1
            if self.logLevel >= LogLevel.General:
                print("Step:{}".format(self.steps))
            res = self.OneStep()
        return res

    def Integrate(self) -> [AStarResult, complex]:
        res: AStarResult = self.FindPath()
        resV: complex = 0
        if AStarResult.Finished == res:
            endNote = self.gridArray[self.endY][self.endX]
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                resV = resV + endNote.neighbourV[endNote.parentDir]
                endNote = endNote.parent
            resV = -resV
            return [res, resV]
        return [AStarResult.Failed, cmath.nan]

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
