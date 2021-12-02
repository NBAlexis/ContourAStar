import math

from UsefulFunctions.NestedQuadrature import NestedQuadrature


class ClenshawCurtisExp(NestedQuadrature):

    def __init__(self, smallProtect: float = 1.0e-15):
        super().__init__()
        self.smallProtect = smallProtect

    def CachePointAndWeightsOneOrder(self, n: int) -> [list, list]:
        """
        see:
        https://people.sc.fsu.edu/~jburkardt/cpp_src/clenshaw_curtis_rule/clenshaw_curtis_rule.html
        """
        if 1 == n:
            return [[0.0], [2.0]]
        totalInterval: int = (1 << (n - 1))
        totalNumber: int = totalInterval + 1
        sep = 1.0 / totalInterval
        theta = [math.pi * (1.0 - sep * i) for i in range(0, totalNumber)]
        w = [1.0 for i in range(0, totalNumber)]
        upper = totalInterval // 2
        for i in range(0, totalNumber):
            t = i * math.pi / totalInterval
            for j in range(1, upper + 1):
                b = 2.0
                if (2 * j) == (totalNumber - 1):
                    b = 1.0
                w[i] = w[i] - b * math.cos(2.0 * j * t) / (4 * j * j - 1)
        w[0] = w[0] / totalInterval
        for i in range(1, totalInterval):
            w[i] = 2.0 * w[i] / totalInterval
        w[totalInterval] = w[totalInterval] / totalInterval
        p = [math.cos(theta[i]) for i in range(0, totalNumber)]
        p[0] = -1.0 + self.smallProtect
        p[totalInterval] = 1.0 - self.smallProtect
        p[totalInterval // 2] = 0.0
        return [p, w]

    def GetNewPoint(self, n: int) -> list:
        if 1 == n:
            return [0]
        if 2 == n:
            return [-1 + self.smallProtect, 1 - self.smallProtect]
        pointNumber: int = 1 << (n - 2)
        sep = 1.0 / (1 << (n - 3))
        start = sep / 2.0
        p = []
        for i in range(0, pointNumber):
            p.append(math.cos(math.pi * (1.0 - start - i * sep)))
        return p

    def GetPointWeightIndex(self, n: int, batch: int, index: int) -> int:
        if batch > n:
            return -1
        totalNumber: int = (1 << (n - 1)) + 1
        if 1 == batch:
            if 1 == n:
                return 0
            return totalNumber // 2
        if 2 == batch:
            return 0 if (0 == index) else (totalNumber - 1)
        sep: int = 1 << (n - batch + 1)
        start = sep // 2
        return start + index * sep

    def GetLeftMostPoint(self) -> [float, int, int]:
        return [-1 + self.smallProtect, 2, 0]


def TestClenshawCurtisExpWeightList(maxOrder: int = 8):
    quadrature = ClenshawCurtisExp()
    print(quadrature.GetPointAndWeight(2))
    print(quadrature.GetPointAndWeight(4))
    print(quadrature.GetPointAndWeight(5))
    [points, starts, ends, weightList] = quadrature.ConstructPointListYLeftMost(maxOrder)
    print(points[0].x)
    print(starts[0])
    print(starts[1])
    print(ends[0])
    print(ends[1])
    print(weightList[0])
    print(weightList[1])
    print(weightList[2])
    print(weightList[3])
