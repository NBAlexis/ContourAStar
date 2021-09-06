"""
10.1023/A:1019129717644
我们暂时只做3D的Gauss-Patterson

第一个点：1,1,1
从2开始：
新增：1,1,2 + 1,2,1 + 2,1,1
新增：1,1,3 + 1,3,1 + 3,1,1 + 1,2,2 + 2,2,1 + 2,1,2
新增：1,1,4 + 1,4,1 + 4,1,1 + 1,2,3 + 1,3,2 + 2,3,1 + 3,2,1 + 2,1,3 + 3,1,2 + 2,2,2

这是：
for i in range(1, n):
 for j in range(1, n - i + 1):
  (i,j,n-i-j+1)

n = 2:
(1,1,1)

n = 3:
i = 1: j in range(1, 3) 1,1,2, 1,2,1
i = 2: j in range(1, 2) 2,1,1

n = 4:
i = 1: j in range(1, 4) 1,1,3, 1,2,2, 1,3,1
i = 2: j in range(1, 3) 2,1,2, 2,2,1
i = 3: j in range(1, 2) 3,1,1

"""
from UsefulFunctions.GaussianPatterson import *
from UsefulFunctions.SparseGridPoints import SparseGridPoints3D


class SparseGrid3D:
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

    def GetWeight(self,
                  xBatch: int, xIndex: int,
                  yBatch: int, yIndex: int,
                  zBatch: int, zIndex: int, n: int) -> float:
        """
        假定 wi是Qi中，这个点的weight
        如果我们规定，点在Qi里新增，则w(j<i)=0
        那么可以定义 di = wi - w(i-1)
        n = 1: d1x * d1y * d1z ...
        n = 2: d1x * d1y * d1z + d2x * d1y * d1z + d1x * d2y * d1z + d1x * d1y * d2z
        ...
        所以可以看做是从m in range(1, n + 1)
        dmx * dmy * dmz的求和。
        dmx按照加点时的batch来。
        """
        res: float = 0
        for m in range(2, n + 2):
            for i in range(1, m):
                for j in range(1, m - i + 1):
                    xList = i
                    yList = j
                    zList = m - i - j + 1
                    dx: float = self.GetDi(xList, xBatch, xIndex)
                    dy: float = self.GetDi(yList, yBatch, yIndex)
                    dz: float = self.GetDi(zList, zBatch, zIndex)
                    res = res + dx * dy * dz
                    # print("adding {}, {}, {}".format(dx, dy, dz))
        return res

    def GetSparseGridNewPoints(self, n: int) -> list:
        """
        假定 n <= 9
        """
        Ql = []
        for i in range(1, n + 1):
            Ql.append(self.quadrature.GetNewPoint(i))
        res = []
        m = n + 1
        for i in range(1, m):
            for j in range(1, m - i + 1):
                xList = i - 1
                yList = j - 1
                zList = m - i - j
                for xi in range(0, len(Ql[xList])):
                    for yi in range(0, len(Ql[yList])):
                        for zi in range(0, len(Ql[zList])):
                            res.append(
                                SparseGridPoints3D(
                                    Ql[xList][xi],
                                    Ql[yList][yi],
                                    Ql[zList][zi],
                                    xList + 1, xi,
                                    yList + 1, yi,
                                    zList + 1, zi))
        return res

    def FillWeights(self, n: int, points: list):
        for point in points:
            point.w = self.GetWeight(
                point.xBatch, point.xIndex,
                point.yBatch, point.yIndex,
                point.zBatch, point.zIndex, n)

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


def TestSparseGridPoints3D(n: int, nestedQuadrature=GaussPatterson()):
    sparseGrid = SparseGrid3D(nestedQuadrature)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(1, n + 1):
        xPoints = []
        yPoints = []
        zPoints = []
        newPoints = sparseGrid.GetSparseGridNewPoints(i)
        for point in newPoints:
            xPoints.append(point.x)
            yPoints.append(point.y)
            zPoints.append(point.z)
        ax.scatter(xPoints, yPoints, zPoints)
    plt.show()


def TestGridWithWeight3D(n: int, nestedQuadrature=GaussPatterson()):
    sparseGrid = SparseGrid3D(nestedQuadrature)
    points = []
    for i in range(1, n + 1):
        points = points + sparseGrid.GetSparseGridNewPoints(i)
    sparseGrid.FillWeights(n, points)
    import matplotlib.pyplot as plt
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for point in points:
        ax.text(point.x, point.y, point.z, "{:.2f}".format(point.w), ha="center")
    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_zlim([-1.2, 1.2])
    plt.show()
