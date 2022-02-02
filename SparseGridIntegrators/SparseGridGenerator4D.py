"""
10.1023/A:1019129717644

第一个点：1,1,1,1
从2开始：
新增：1,1,1,2 + 1,1,2,1 + 1,2,1,1 + 2,1,1,1
新增：1,1,1,3 + 1,1,3,1 + 1,3,1,1 + 3.1,1,1 + 2,1,1,2 + 2,1,2,1 + 2,2,1,1 + 1,2,1,2 + 1,2,2,1 + 1,1,2,2
新增：1,1,1,4 + ... +  + 1,3,2 + 2,3,1 + 3,2,1 + 2,1,3 + 3,1,2 + 2,2,2

这是：
n=3
for m in range(1, n+1):
    for i in range(1, m+1):
        for j in range(1, m - i + 2):
            print(i , j, m-i-j+2, n-m+1)

n = 1:
1 1 1 1

n = 2:
1 1 1 2
1 1 2 1
1 2 1 1
2 1 1 1

n = 3:
1 1 1 3
1 1 2 2
1 2 1 2
2 1 1 2
1 1 3 1
1 2 2 1
1 3 1 1
2 1 2 1
2 2 1 1
3 1 1 1

"""
from SparseGridIntegrators.NestedQuadrature import NestedQuadrature
from SparseGridIntegrators.SparseGridPoints import SparseGridPoints4D


class SparseGrid4D:
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
                  zBatch: int, zIndex: int,
                  wBatch: int, wIndex: int, n: int) -> float:
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
        for nn in range(1, n + 1):
            for m in range(1, nn + 1):
                for i in range(1, m + 1):
                    for j in range(1, m - i + 2):
                        # print(i, j, m - i - j + 2, nn - m + 1)
                        xList = i
                        yList = j
                        zList = m - i - j + 2
                        wList = nn - m + 1
                        dx: float = self.GetDi(xList, xBatch, xIndex)
                        dy: float = self.GetDi(yList, yBatch, yIndex)
                        dz: float = self.GetDi(zList, zBatch, zIndex)
                        dw: float = self.GetDi(wList, wBatch, wIndex)
                        res = res + dx * dy * dz * dw
                        # print("adding {}, {}, {}, {}".format(dx, dy, dz, dw))
        return res

    def GetSparseGridNewPoints(self, n: int) -> list:
        """
        假定 n <= 9
        """
        Ql = []
        for i in range(1, n + 1):
            Ql.append(self.quadrature.GetNewPoint(i))
        res = []
        for m in range(1, n + 1):
            for i in range(1, m + 1):
                for j in range(1, m - i + 2):
                    xList = i - 1
                    yList = j - 1
                    zList = m - i - j + 1
                    wList = n - m
                    for xi in range(0, len(Ql[xList])):
                        for yi in range(0, len(Ql[yList])):
                            for zi in range(0, len(Ql[zList])):
                                for wi in range(0, len(Ql[wList])):
                                    res.append(
                                        SparseGridPoints4D(
                                            Ql[xList][xi],
                                            Ql[yList][yi],
                                            Ql[zList][zi],
                                            Ql[wList][wi],
                                            xList + 1, xi,
                                            yList + 1, yi,
                                            zList + 1, zi,
                                            wList + 1, wi))
        return res

    def FillWeights(self, n: int, points: list):
        for point in points:
            point.weight = self.GetWeight(
                point.xBatch, point.xIndex,
                point.yBatch, point.yIndex,
                point.zBatch, point.zIndex,
                point.wBatch, point.wIndex, n)

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
            print("Initial sparse grid integrator: generate order ", n)
            newPoint = self.GetSparseGridNewPoints(n)
            starts.append(len(points))
            points = points + newPoint
            ends.append(len(points))
            self.FillWeights(n, points)
            weights = []
            for point in points:
                weights.append(point.weight)
            weightList.append(weights)
        return [points, starts, ends, weightList]
