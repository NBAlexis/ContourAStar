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
