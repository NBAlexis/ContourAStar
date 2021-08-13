from UsefulFunctions.SparseGridPoints import SparseGridPoints


class NestedQuadrature:
    def GetPointAndWeight(self, n: int) -> [list, list]:
        pass

    def GetNewPoint(self, n: int) -> list:
        pass

    def GetPointWeightIndex(self, n: int, batch: int, index: int) -> int:
        pass

    def ConstructPointListY(self, xBatch: int, xIndex: int, maxOrder: int = 12) -> [list, list, list, list]:
        points = []
        starts = []
        ends = []
        weightList = []
        [px, _] = self.GetPointAndWeight(xBatch)
        xCoord = px[((xIndex + 1) << 1) - 2]
        if self.MaxOrder() > 0:
            maxOrder = min(maxOrder, self.MaxOrder())
        for n in range(1, maxOrder + 1):
            yCoords = self.GetNewPoint(n)
            newPoint = []
            for yIndex in range(0, len(yCoords)):
                newPoint.append(SparseGridPoints(xCoord, yCoords[yIndex], xBatch, xIndex, n, yIndex))
            starts.append(len(points))
            points = points + newPoint
            ends.append(len(points))
            newWeights = []
            [_, wy] = self.GetPointAndWeight(n)
            for point in points:
                newWeights.append(wy[self.GetPointWeightIndex(n, point.yBatch, point.yIndex)])
            weightList.append(newWeights.copy())
        return [points, starts, ends, weightList]

    def ConstructPointListYLeftMost(self, maxOrder: int) -> [list, list, list, list]:
        [_, batchX, indexX] = self.GetLeftMostPoint()
        return self.ConstructPointListY(batchX, indexX, maxOrder)

    def GetLeftMostPoint(self) -> [float, int, int]:
        pass

    def MaxOrder(self) -> int:
        return -1


class Trapezoidal(NestedQuadrature):

    def __init__(self, smallProtect: float = 1.0e-15):
        self.smallProtect = smallProtect

    def GetPointAndWeight(self, n: int) -> [list, list]:
        pointNumber: int = 1 << (n + 1)
        sep: float = 2 / pointNumber
        p = []
        w = []
        for i in range(0, pointNumber + 1):
            thisP = i * sep - 1
            thisW = sep / 3
            if 0 == i:
                thisP = -1 + self.smallProtect
            elif pointNumber == i:
                thisP = 1 - self.smallProtect
            elif 0 == (i & 1):
                thisW = 2 * sep / 3
            else:
                thisW = 4 * sep / 3
            p.append(thisP)
            w.append(thisW)
        return [p, w]

    def GetNewPoint(self, n: int) -> list:
        if 1 == n:
            return [-1 + self.smallProtect, -0.5, 0, 0.5, 1 - self.smallProtect]
        pointNumber: int = 1 << n
        start: float = 1 / pointNumber
        sep: float = 2 * start
        p = []
        for i in range(0, pointNumber):
            p.append(start + i * sep - 1)
        return p

    def GetPointWeightIndex(self, n: int, batch: int, index: int) -> int:
        if batch > n:
            return -1
        if 1 == batch:
            return index << (n - 1)
        # just Gauss Patterson add one
        return (((index + 1) << 1) - 1) << (n - batch)

    def GetLeftMostPoint(self) -> [float, int, int]:
        return [-1 + self.smallProtect, 1, 0]


def TestTrapezoidalWeightList(maxOrder: int = 12):
    quadrature = Trapezoidal()
    print(quadrature.GetPointAndWeight(3))
    [points, starts, ends, weightList] = quadrature.ConstructPointListYLeftMost(maxOrder)
    print(points[0].x)
    print(starts[0])
    print(starts[1])
    print(ends[0])
    print(ends[1])
    print(weightList[0])
    print(weightList[1])
    print(weightList[2])
