"""
10.1023/A:1019129717644
我们暂时只做2D的Gauss-Patterson

第一个点：1,1
从2开始：
新增：1,2 + 2,1
新增：1,3 + 2,2 + 3,1
新增：1,4 + 2,3 + 3,2 + 4,1
新增：1,5 + 2,4 + 3,3 + 4,2 + 5,1
新增：1,6 + 2,5 + 3,4 + 4,3 + 5,2 + 6,1
新增：1,7 + 2,6 + 3,5 + 4,4 + 5,3 + 6,2 + 7,1
新增：1,8 + 2,7 + 3,6 + 4,5 + 5,4 + 6,3 + 7,2 + 8,1
新增：1,9 + 2,8 + 3,7 + 4,6 + 5,5 + 6,4 + 7,3 + 8,2 + 9,1

"""
from UsefulFunctions.GaussianPatterson import *
from UsefulFunctions.SparseGridPoints import SparseGridPoints


class SparseGrid:
    def __init__(self, quadrature: NestedQuadrature):
        self.quadrature = quadrature

    def GetDi(self, batch: int, xBatch: int, xIndex: int) -> float:
        """
        注意！这个函数只在构造weight的时候调用，因为它很卡
        """
        weightI: float = 0
        weightIm1: float = 0
        wi = self.quadrature.GetPointWeightIndex(batch, xBatch, xIndex)
        wim1 = self.quadrature.GetPointWeightIndex(batch - 1, xBatch, xIndex)
        if wi >= 0:
            [_, w1] = self.quadrature.GetPointAndWeight(batch)
            weightI = w1[wi]
        if wim1 >= 0:
            [_, w2] = self.quadrature.GetPointAndWeight(batch - 1)
            weightIm1 = w2[wim1]
        return weightI - weightIm1

    def GetWeight(self, xBatch: int, xIndex: int, yBatch: int, yIndex: int, n: int) -> float:
        """
        假定 wi是Qi中，这个点的weight
        如果我们规定，点在Qi里新增，则w(j<i)=0
        那么可以定义 di = wi - w(i-1)
        n = 1: d1x * d1y
        n = 2: d1x * d1y + d2x * d1y + d1x * d2y
        n = 3: d1x * d1y + d2x * d1y + d1x * d2y + d3x * d1y + d2x * d2y + d1x * d3y
        ...
        """
        res: float = 0
        for i in range(1, n + 1):
            for j in range(1, n + 2 - i):
                dix: float = self.GetDi(i, xBatch, xIndex)
                djy: float = self.GetDi(j, yBatch, yIndex)
                res = res + dix * djy
        return res

    def GetSparseGridNewPoints(self, n: int) -> list:
        """
        假定 n <= 9
        """
        Ql = []
        for i in range(1, n + 1):
            Ql.append(self.quadrature.GetNewPoint(i))
        res = []
        for i in range(1, n + 1):
            xList = i - 1
            yList = n - i
            for xi in range(0, len(Ql[xList])):
                for yi in range(0, len(Ql[yList])):
                    res.append(SparseGridPoints(Ql[xList][xi], Ql[yList][yi], xList + 1, xi, yList + 1, yi))
        return res

    def FillWeights(self, n: int, points: list):
        for point in points:
            point.w = self.GetWeight(point.xBatch, point.xIndex, point.yBatch, point.yIndex, n)

    def ConstructPointListAndWeightList(self, order: int = 20) -> [list, list, list, list]:
        """
        Gauss Patterson 一直创建到order9
        """
        maxOrder = self.quadrature.MaxOrder()
        order = min(order, maxOrder) if maxOrder > 0 else order
        points = []
        starts = []
        ends = []
        weightList = []
        for n in range(1, order + 1):
            print("generate order:", n)
            newPoint = self.GetSparseGridNewPoints(n)
            starts.append(len(points))
            points = points + newPoint
            ends.append(len(points))
            self.FillWeights(n, points)
            weights = []
            for point in points:
                weights.append(point.w)
            weightList.append(weights)
        return [points, starts, ends, weightList]


def TestSparseGridPoints(n: int, nestedQuadrature=GaussPatterson()):
    sparseGrid = SparseGrid(nestedQuadrature)
    import matplotlib.pyplot as plt
    for i in range(1, n + 1):
        xPoints = []
        yPoints = []
        newPoints = sparseGrid.GetSparseGridNewPoints(i)
        for point in newPoints:
            xPoints.append(point.x)
            yPoints.append(point.y)
        plt.scatter(xPoints, yPoints)
    plt.show()


def TestGridWithWeight(n: int, nestedQuadrature=GaussPatterson()):
    sparseGrid = SparseGrid(nestedQuadrature)
    points = []
    for i in range(1, n + 1):
        points = points + sparseGrid.GetSparseGridNewPoints(i)
    sparseGrid.FillWeights(n, points)
    import matplotlib.pyplot as plt
    for point in points:
        plt.text(point.x, point.y, "{:.3f}".format(point.w), ha="center")
    plt.xlim([-1.2, 1.2])
    plt.ylim([-1.2, 1.2])
    plt.show()
