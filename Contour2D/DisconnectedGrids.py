import cmath

from Contour1D.CGrids import GridNeighbour, GridState, AStarResult
from Contour1D.CommonDefinitions import LogLevel
from Contour2D import Integrator2D
from Contour2D.Integrand2D import Integrand2D

"""
This is a proof of concept try
It doesn't work...
"""

class OneDiscGrid:

    def __init__(self, owner, coords: list, values: list, idx: int, neighbourCount: int):
        self.owner = owner
        self.parent = None
        self.parentDir = -1
        self.coords = coords
        self.values = values
        self.idx = idx
        self.neighbourV = [0j for i in range(0, neighbourCount)]
        self.neighbourNode = [None for i in range(0, neighbourCount)]
        self.neighbourState = [GridNeighbour.Unknown for i in range(0, neighbourCount)]
        self.state = GridState.NotDecide
        self.fn = 0
        self.hn = 0
        self.gn = 0
        self.start = False
        self.target = False

    def SetNeighbour(self, direction: int, neighbour):
        if neighbour is None:
            self.neighbourState[direction] = GridNeighbour.NotConnected
            return
        self.neighbourNode[direction] = neighbour

    """
    def UpdateNeighbour(self, direction: int, integrateValue: complex, state: GridNeighbour):
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
    """

    def SetHn(self, targetX, targetY):
        self.hn = abs(self.coords[0] - targetX) + abs(self.coords[1] - targetY) + abs(self.coords[2] - targetX) + abs(self.coords[3] - targetY)

    def UpdateFn(self, parentNode, direction: int):
        oppositeD: int = direction + 1 - ((direction & 1) << 1)
        if self.state == GridState.NotDecide:
            self.parent = parentNode
            self.parentDir = oppositeD
            self.gn = 1 + parentNode.gn
            self.fn = self.hn + self.gn
            self.state = GridState.OpenList
            return
        newGn = 1 + parentNode.gn
        if newGn < self.gn:
            self.gn = newGn
            self.parent = parentNode
            self.parentDir = oppositeD
            self.fn = self.gn + self.hn


class DisconnectedGrids2D:
    """
    假设width, height都是大于等于3的整数
    假设edge < (width - 1) / 2
    """
    def __init__(self, width: int, edge: int, integrator: Integrator2D, integrand: Integrand2D,
                 maxStep: int = 1000, logLevel: LogLevel = LogLevel.General):
        self.integrator = integrator
        self.integrand = integrand
        midY = width // 2
        self.sep = 2.0 / (width - 1 - 2 * edge)
        self.width = width
        self.edge = edge
        self.startX = edge
        self.startY = midY
        self.endX = width - 1 - edge
        self.endY = midY
        self.maxStep = maxStep
        self.logLevel = logLevel
        self.steps = 0
        self.AStarRes = AStarResult.Unknown
        self.gridArray = [[[[
            OneDiscGrid(self, [rex, imx, rey, imy],
                        [-1 + (rex - edge) * self.sep + (imx - midY) * self.sep * 1j,
                         -1 + (rey - edge) * self.sep + (imy - midY) * self.sep * 1j],
                        rex + imx * width + rey * width * width + imy * width * width * width,
                        8)
            for rex in range(0, width)]
            for imx in range(0, width)]
            for rey in range(0, width)]
            for imy in range(0, width)]
        self.gridList = []
        """
        方向： 0,1,2,3,4,5,6,7 分别是实轴x+,实轴x-，虚轴x+,...
        反方向：偶数+1， 奇数-1
        也就是opposite = dir + 1 - ((dir & 1) << 1)
        注意gridArray是[imy][rey][imx][rex]
        offset: 
        """
        rexOffsets = [1, -1, 0, 0, 0, 0, 0, 0]
        imxOffsets = [0, 0, 1, -1, 0, 0, 0, 0]
        reyOffsets = [0, 0, 0, 0, 1, -1, 0, 0]
        imyOffsets = [0, 0, 0, 0, 0, 0, 1, -1]
        for imy in range(0, width):
            for rey in range(0, width):
                for imx in range(0, width):
                    for rex in range(0, width):
                        self.gridArray[imy][rey][imx][rex].SetHn(self.endX, self.endY)
                        self.gridList.append(self.gridArray[imy][rey][imx][rex])
                        # fix up the neighbours
                        for d in range(0, 8):
                            rexOffset: int = rexOffsets[d] + rex
                            imxOffset: int = imxOffsets[d] + imx
                            reyOffset: int = reyOffsets[d] + rey
                            imyOffset: int = imyOffsets[d] + imy
                            if 0 <= rexOffset < self.width and 0 <= imxOffset < self.width and 0 <= reyOffset < self.width and 0 <= imyOffset < self.width:
                                self.gridArray[imy][rey][imx][rex].SetNeighbour(d, self.gridArray[imyOffset][reyOffset][imxOffset][rexOffset])
                            else:
                                self.gridArray[imy][rey][imx][rex].SetNeighbour(d, None)
        self.ResetGrid()

    def ResetGrid(self):
        for imy in range(0, self.width):
            for rey in range(0, self.width):
                for imx in range(0, self.width):
                    for rex in range(0, self.width):
                        self.gridArray[imy][rey][imx][rex].SetHn(self.endX, self.endY)
                        self.gridArray[imy][rey][imx][rex].parent = None
                        self.gridArray[imy][rey][imx][rex].parentDir = -1
                        self.gridArray[imy][rey][imx][rex].fn = 0
                        self.gridArray[imy][rey][imx][rex].state = GridState.NotDecide
                        for d in range(0, 8):
                            if self.gridArray[imy][rey][imx][rex].neighbourNode[d] is None:
                                self.gridArray[imy][rey][imx][rex].neighbourState[d] = GridNeighbour.NotConnected
                            else:
                                self.gridArray[imy][rey][imx][rex].neighbourState[d] = GridNeighbour.Unknown
        self.gridArray[self.startY][self.startX][self.startY][self.startX].state = GridState.OpenList
        self.gridArray[self.startY][self.startX][self.startY][self.startX].start = True
        self.gridArray[self.endY][self.endX][self.endY][self.endX].target = True
        self.steps = 0
        self.AStarRes = AStarResult.Unknown

    def CalculateConnection(self, node: OneDiscGrid, direction: int):
        if node.neighbourNode[direction] is None:
            node.neighbourState[direction] = GridNeighbour.NotConnected
            return
        if node.neighbourState[direction] == GridNeighbour.Connected:
            return
        oppositeD: int = direction + 1 - ((direction & 1) << 1)
        bHasValue = False
        value = 0j
        if direction <= 3:
            # y is not changed, only change x
            print("x: ", node.values[0], " ", node.neighbourNode[direction].values[0], "y: ", node.values[1], " 1")
            [bHasValue, value] = self.integrator.Integrate(
                self.integrand,
                node.values[0],
                node.neighbourNode[direction].values[0],
                node.values[1], 1)
        else:
            # x is not changed, only change y
            print("x: ", node.values[0], " 1, y: ", node.values[1], " ", node.neighbourNode[direction].values[1])
            [bHasValue, value] = self.integrator.Integrate(
                self.integrand,
                node.values[0], 1,
                node.values[1],
                node.neighbourNode[direction].values[1])
        if bHasValue:
            node.neighbourState[direction] = GridNeighbour.Connected
            node.neighbourV[direction] = value
            node.neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.Connected
            node.neighbourNode[direction].neighbourV[oppositeD] = -value
        else:
            node.neighbourState[direction] = GridNeighbour.NotConnected
            node.neighbourNode[direction].neighbourState[oppositeD] = GridNeighbour.NotConnected

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
        for i in range(0, 8):
            if smallestFnNode.neighbourState[i] == GridNeighbour.Unknown:
                self.CalculateConnection(smallestFnNode, i)
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
            endNote = self.gridArray[self.endY][self.endX][self.endY][self.endX]
            while (endNote is not None) and (not endNote.start):
                if self.logLevel >= LogLevel.General:
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
                for d in range(0, 8):
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