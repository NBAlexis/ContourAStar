import math

from SparseGridIntegrators.NestedQuadrature import NestedQuadrature


class ClenshawCurtis(NestedQuadrature):
    """
    https://people.math.sc.edu/Burkardt/f77_src/ccn_rule/ccn_rule.html

    """
    def __init__(self, smallProtect: float = 1.0e-15):
        super().__init__()
        self.smallProtect = smallProtect

    def GetNewPoint(self, n: int) -> list:
        """
        1/2,
        0, 1
        1/4, 3/4
        1/8, 3/8, 5/8, 7/8,
        1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
        But we would prefer that the numbers in each row be regrouped in pairs
        that are symmetric about 1/2, with the number above 1/2 coming first.
        Thus, the last row might become:
        (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
        """
        if 1 == n:
            return [0]
        if 2 == n:
            return [-1 + self.smallProtect, 1 - self.smallProtect]
        # 2, 4, 8 pairs ...
        pair = -1
        pairIndex = 0
        orderCount = n - 2
        while orderCount > 0:
            pair = pair + 1
            pairCount = 1 << pair
            pairIndex = 0
            for i in range(0, pairCount):
                orderCount = orderCount - 1
                pairIndex = pairIndex + 1
                if 0 == orderCount:
                    break
        # for n = 3, it will stop at pair:0, pairIndex:1
        #         4                       1            1
        #         5                       1            2
        #         6                       2            1
        #         7                       2            2
        #         8                       2            3
        #         9                       2            4
        #        10                       3            1
        mid = 1 << (pair + 1)
        denominator = (mid << 1)
        left = mid + (2 * pairIndex - 1)
        right = mid - (2 * pairIndex - 1)
        # for n = 3, it will stop at pair:0, pairIndex:1   mid: 2  denominator: 4  left: 3  right: 1
        #         4                       1            1        4               8        5         3
        #         5                       1            2        4               8        7         1
        #         6                       2            1        8              16        9         7
        #         7                       2            2        8              16       11         5
        #         8                       2            3        8              16
        #         9                       2            4        8              16
        #        10                       3            1
        return [math.cos(math.pi * left / denominator), math.cos(math.pi * right / denominator)]

    def GetPointWeightIndex(self, n: int, batch: int, index: int) -> int:
        if batch > n:
            return -1
        if 1 == batch:
            return 0
        return (batch << 1) - 3 + index

    def GetLeftMostPoint(self) -> [float, int, int]:
        return [-1 + self.smallProtect, 2, 0]

    def CachePointAndWeightsOneOrder(self, toCache: int):
        """
        1/2,
        0, 1
        1/4, 3/4
        1/8, 3/8, 5/8, 7/8,
        1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16, and so on.
        But we would prefer that the numbers in each row be regrouped in pairs
        that are symmetric about 1/2, with the number above 1/2 coming first.
        Thus, the last row might become:
        (9/16, 7/16), (11/16, 5/16), ..., (15/16, 1/16).
        """
        order = (toCache << 1) - 1
        p = []
        for i in range(1, toCache + 1):
            p = p + self.GetNewPoint(i)
        w = [0 for _ in range(0, order)]
        for i in range(0, order):
            # Compute the Lagrange basis polynomial which is 1 at XTAB(I), and zero at the other nodes.
            d = [0 for _ in range(0, order)]
            d[i] = 1
            for j in range(2, order + 1):
                for k in range(j, order + 1):
                    idx = order + j - k - 1
                    d[idx] = (d[idx - 1] - d[idx]) / (p[order - k] - p[idx])
            for j in range(1, order):
                for k in range(1, order - j + 1):
                    idx = order - k - 1
                    d[idx] = d[idx] - p[order - k - j] * d[idx + 1]
            # Evaluate the antiderivative of the polynomial at the left and right endpoints.
            ya = d[order - 1] / order
            for j in range(order - 2, -1, -1):
                ya = -ya + d[j] / (j + 1)
            ya = -ya
            yb = d[order - 1] / order
            for j in range(order - 2, -1, -1):
                yb = yb + d[j] / (j + 1)
            w[i] = yb - ya
        return [p, w]


def TestClenshawCurtisWeightList(maxOrder: int = 12):
    quadrature = ClenshawCurtis()
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

