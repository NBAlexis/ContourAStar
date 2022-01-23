import cmath
import random

from Contour1D.CGrids import GridDir, GridNeighbour, GridState, OneGrid, AStarResult
from Contour1D.CommonDefinitions import LogLevel


def TraceBack(checkGrid: OneGrid) -> list:
    retList = []
    startNode = checkGrid
    parentNode = checkGrid.parent
    while parentNode is not None:
        retList.append([parentNode, startNode])
        startNode = parentNode
        parentNode = parentNode.parent
    return retList


class CPath:

    def __init__(self, nodes: list, points: list = None):
        self.nodes = nodes
        self.points = points

    def EqualTo(self, other) -> bool:
        if len(self.nodes) != len(other.nodes):
            return False
        for i in range(0, len(self.nodes)):
            if self.nodes[i].idx != other.nodes[i].idx:
                return False
        return True


class CGeneralGrid:

    def __init__(self, width: int, height: int, edge: int,
                 maxStep: int = 1000000, logLevel: LogLevel = LogLevel.General):
        self.midY = height // 2
        self.sep = 2.0 / (width - 1 - 2 * edge)
        self.width = width
        self.height = height
        self.edge = edge
        self.startX = edge
        self.startY = self.midY
        self.endX = width - 1 - edge
        self.endY = self.midY
        [self.gridArray, self.gridList] = self.CreateGrid()
        self.steps = 0
        self.maxStep = maxStep
        self.logLevel = logLevel
        self.AStarRes = AStarResult.Unknown
        self.lastError = ""

    def CreateGrid(self) -> [list, list]:
        retList = []
        retGrid = [
            [OneGrid(self, x, y, x + y * self.width, -1 + (x - self.edge) * self.sep + (y - self.midY) * self.sep * 1j)
             for x in range(0, self.width)] for y
            in range(0, self.height)]
        xOffsets = [0, -1, 1, 0]
        yOffsets = [1, 0, 0, -1]
        for y in range(0, self.height):
            for x in range(0, self.width):
                retGrid[y][x].SetHn(self.endX, self.endY)
                retList.append(retGrid[y][x])
                # fix up the neighbours
                for d in range(0, GridDir.Max):
                    xOffset: int = xOffsets[d] + x
                    yOffset: int = yOffsets[d] + y
                    if 0 <= xOffset < self.width and 0 <= yOffset < self.height:
                        retGrid[y][x].SetNeighbour(GridDir(d), retGrid[yOffset][xOffset])
                    else:
                        retGrid[y][x].SetNeighbour(GridDir(d), None)
        retGrid[self.startY][self.startX].start = True
        retGrid[self.endY][self.endX].target = True
        return [retGrid, retList]

    """
    完全重置
    """

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

    """
    不改变连接状态，重置
    """

    def RestartGrid(self):
        for y in range(0, self.height):
            for x in range(0, self.width):
                self.gridArray[y][x].SetHn(self.endX, self.endY)
                self.gridArray[y][x].parent = None
                self.gridArray[y][x].parentDir = GridDir.Max
                self.gridArray[y][x].fn = 0
                self.gridArray[y][x].state = GridState.NotDecide
        self.gridArray[self.startY][self.startX].state = GridState.OpenList
        self.gridArray[self.startY][self.startX].start = True
        self.gridArray[self.endY][self.endX].target = True
        self.steps = 0
        self.AStarRes = AStarResult.Unknown

    def TurnConnectedToUnknown(self):
        for y in range(0, self.height):
            for x in range(0, self.width):
                for i in range(0, len(self.gridArray[y][x].neighbourState)):
                    if self.gridArray[y][x].neighbourState[i] == GridNeighbour.Connected:
                        self.gridArray[y][x].neighbourState[i] = GridNeighbour.Unknown

    def CalculateConnectionStep(self, checkConnectionFunction, direction: GridDir, x: int, y: int):
        if self.gridArray[y][x].neighbourNode[direction] is None:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            return
        if self.gridArray[y][x].neighbourState[direction] == GridNeighbour.Connected:
            return
        if self.gridArray[y][x].neighbourState[direction] == GridNeighbour.NotConnected:
            return
        oppositeD: GridDir = GridDir.Down - direction
        [bHasValue, v] = checkConnectionFunction(self.gridArray[y][x], self.gridArray[y][x].neighbourNode[direction])
        if bHasValue:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourV[direction] = v
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            self.gridArray[y][x].neighbourNode[direction].neighbourV[oppositeD] = -v
        else:
            self.gridArray[y][x].neighbourState[direction] = GridNeighbour.NotConnected
            self.gridArray[y][x].neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

    def GetBlockedInterval(self, alreadyBlocked: list) -> list:
        randomList = []
        for oneGrid in self.gridList:
            if oneGrid.state == GridState.ClosedList:
                for i in range(0, len(oneGrid.neighbourState)):
                    if oneGrid.neighbourState[i] == GridNeighbour.NotConnected and oneGrid.neighbourNode[i] is not None:
                        v1 = oneGrid.v
                        v2 = oneGrid.neighbourNode[i].v
                        hasDuplicated = False
                        for already in alreadyBlocked:
                            if abs(already[0] - v1) + abs(already[1] - v2) < 0.000000001:
                                hasDuplicated = True
                                break
                            if abs(already[0] - v2) + abs(already[1] - v1) < 0.000000001:
                                hasDuplicated = True
                                break
                        if not hasDuplicated:
                            randomList.append([oneGrid.v, oneGrid.neighbourNode[i].v])
        if len(randomList) > 0:
            return randomList[random.randint(0, len(randomList) - 1)]
        return []

    def GetBlockedIntervalNodes(self, alreadyBlockedNodes: list) -> list:
        randomList = []
        for oneGrid in self.gridList:
            if oneGrid.state == GridState.ClosedList:
                for i in range(0, len(oneGrid.neighbourState)):
                    if oneGrid.neighbourState[i] == GridNeighbour.NotConnected and oneGrid.neighbourNode[i] is not None:
                        v1: OneGrid = oneGrid
                        v2: OneGrid = oneGrid.neighbourNode[i]
                        hasDuplicated = False
                        for already in alreadyBlockedNodes:
                            if already[0].idx == v1.idx and already[1].idx == v2.idx:
                                hasDuplicated = True
                                break
                            if already[1].idx == v1.idx and already[0].idx == v2.idx:
                                hasDuplicated = True
                                break
                        if not hasDuplicated:
                            randomList.append([v1, v2])
        if len(randomList) > 0:
            return randomList[random.randint(0, len(randomList) - 1)]
        return []

    def ExtendNearestPath(self) -> list:
        nearestDist = 1000000.0
        allFarNodes = []
        allPath = []
        for oneGrid in self.gridList:
            if oneGrid.state == GridState.ClosedList:
                if abs(oneGrid.hn - nearestDist) < 0.0001:
                    allFarNodes.append(oneGrid)
                elif oneGrid.hn < nearestDist:
                    nearestDist = oneGrid.hn
                    allFarNodes = [oneGrid]
        for oneGrid in allFarNodes:
            thisGridPath = TraceBack(oneGrid)
            for i in range(0, len(oneGrid.neighbourState)):
                if oneGrid.neighbourState[i] == GridNeighbour.NotConnected and oneGrid.neighbourNode[i] is not None:
                    toAdd = [thisGridPath[n][0] for n in range(0, len(thisGridPath))]
                    toAdd.append(oneGrid)
                    toAdd.append(oneGrid.neighbourNode[i])
                    allPath.append(CPath(toAdd))
        return allPath

    def ExtendAllPath(self) -> list:
        allPath = []
        for oneGrid in self.gridList:
            if oneGrid.state == GridState.ClosedList:
                thisGridPath = TraceBack(oneGrid)
                for i in range(0, len(oneGrid.neighbourState)):
                    if oneGrid.neighbourState[i] == GridNeighbour.NotConnected and oneGrid.neighbourNode[i] is not None:
                        toAdd = [thisGridPath[n][0] for n in range(0, len(thisGridPath))]
                        toAdd.append(oneGrid)
                        toAdd.append(oneGrid.neighbourNode[i])
                        allPath.append(CPath(toAdd))
        return allPath

    def OneStep(self, checkConnectionFunction) -> AStarResult:
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
                self.CalculateConnectionStep(checkConnectionFunction, GridDir(i), smallestFnNode.x, smallestFnNode.y)
            if smallestFnNode.neighbourState[i] == GridNeighbour.Connected:
                neighbour = smallestFnNode.neighbourNode[i]
                if neighbour.state != GridState.ClosedList:
                    neighbour.UpdateFn(smallestFnNode, i)
        smallestFnNode.state = GridState.ClosedList
        return AStarResult.Unknown

    def FindPath(self, checkConnectionFunction) -> [AStarResult, list, complex]:
        self.RestartGrid()
        res: AStarResult = self.OneStep(checkConnectionFunction)
        while AStarResult.Unknown == res and self.steps < self.maxStep:
            self.steps = self.steps + 1
            if self.logLevel >= LogLevel.Verbose:
                print("Step:{}".format(self.steps))
            res = self.OneStep(checkConnectionFunction)
        resV: complex = 0
        if AStarResult.Finished == res:
            endNote = self.gridArray[self.endY][self.endX]
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                resV = resV + endNote.neighbourV[endNote.parentDir]
                endNote = endNote.parent
            # gather path
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
            return [res, points, -resV]
        return [AStarResult.Failed, [], cmath.nan]

    def FindPathWithNodes(self, checkConnectionFunction) -> [AStarResult, list, list, complex]:
        self.RestartGrid()
        res: AStarResult = self.OneStep(checkConnectionFunction)
        while AStarResult.Unknown == res and self.steps < self.maxStep:
            self.steps = self.steps + 1
            if self.logLevel >= LogLevel.Verbose:
                print("Step:{}".format(self.steps))
            res = self.OneStep(checkConnectionFunction)
        resV: complex = 0
        if AStarResult.Finished == res:
            endNote = self.gridArray[self.endY][self.endX]
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                resV = resV + endNote.neighbourV[endNote.parentDir]
                endNote = endNote.parent
            # gather path
            lastDir = GridDir.Max
            endNote = self.gridArray[self.endY][self.endX]
            points = []
            nodes = []
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.Verbose:
                    print("dir:{}, v:{}".format(endNote.parentDir, endNote.neighbourV[endNote.parentDir]))
                if lastDir != endNote.parentDir:
                    points.insert(0, endNote.v)
                    lastDir = endNote.parentDir
                nodes.insert(0, endNote)
                endNote = endNote.parent
            points.insert(0, endNote.v)
            nodes.insert(0, endNote)
            return [res, nodes, points, -resV]
        return [AStarResult.Failed, [], [], cmath.nan]

    def GetPathNodes(self) -> list:
        endNote = self.gridArray[self.endY][self.endX]
        points = []
        while (endNote is not None) and (not endNote.start):
            points.insert(0, endNote)
            endNote = endNote.parent
        points.insert(0, endNote)
        return points

    def EliminateMiddlePointOfPath(self) -> bool:
        node = []
        endNode = self.gridArray[self.endY][self.endX]
        node.append(endNode)
        while (endNode is not None) and (not endNode.start):
            endNode = endNode.parent
            node.append(endNode)
        if len(node) < 3:
            return False
        toRemove = node[len(node) // 2 - 1]
        toRemove.neighbourState[toRemove.parentDir] = GridNeighbour.NotConnected
        toRemove.neighbourNode[toRemove.parentDir].neighbourState[
            GridDir.Down - toRemove.parentDir] = GridNeighbour.NotConnected

    def Show(self):
        import matplotlib.pyplot as plt
        sep = self.sep * 0.25
        xOffsets = [0, -sep, sep, 0]
        yOffsets = [sep, 0, 0, -sep]
        for y in range(0, self.height):
            for x in range(0, self.width):
                v = self.gridArray[y][x].v
                xp = v.real
                yp = v.imag
                colorList = ['black', 'blue', 'red']
                plt.text(xp, yp, "({},{})".format(x, y), color=colorList[self.gridArray[y][x].state])
                for d in range(0, GridDir.Max):
                    if self.gridArray[y][x].neighbourState[d] == GridNeighbour.Unknown:
                        plt.arrow(x=xp, y=yp, dx=xOffsets[d], dy=yOffsets[d], width=.02, facecolor='yellow')
                    elif self.gridArray[y][x].neighbourState[d] == GridNeighbour.Connected:
                        if self.gridArray[y][x].parentDir == d:
                            plt.arrow(x=xp, y=yp, dx=xOffsets[d], dy=yOffsets[d], width=.02, facecolor='green')
                        else:
                            plt.arrow(x=xp, y=yp, dx=xOffsets[d], dy=yOffsets[d], width=.02, facecolor='blue')
        plt.show()
