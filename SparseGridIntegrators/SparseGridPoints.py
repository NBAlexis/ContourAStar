class SparseGridPoints:
    def __init__(self, x: float, y: float, xBatch: int, xIndex: int, yBatch: int, yIndex: int,
                 w: float = 0, v: complex = 0):
        self.x = x
        self.y = y
        self.w = w
        self.v = v
        self.xBatch = xBatch
        self.xIndex = xIndex
        self.yBatch = yBatch
        self.yIndex = yIndex


class SparseGridPoints3D:
    def __init__(self, x: float, y: float, z: float,
                 xBatch: int, xIndex: int,
                 yBatch: int, yIndex: int,
                 zBatch: int, zIndex: int,
                 w: float = 0, v: complex = 0):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.v = v
        self.xBatch = xBatch
        self.xIndex = xIndex
        self.yBatch = yBatch
        self.yIndex = yIndex
        self.zBatch = zBatch
        self.zIndex = zIndex


class SparseGridPoints4D:
    def __init__(self, x: float, y: float, z: float, w: float,
                 xBatch: int, xIndex: int,
                 yBatch: int, yIndex: int,
                 zBatch: int, zIndex: int,
                 wBatch: int, wIndex: int,
                 weight: float = 0, v: complex = 0):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.weight = w
        self.v = v
        self.xBatch = xBatch
        self.xIndex = xIndex
        self.yBatch = yBatch
        self.yIndex = yIndex
        self.zBatch = zBatch
        self.zIndex = zIndex
        self.wBatch = wBatch
        self.wIndex = wIndex
