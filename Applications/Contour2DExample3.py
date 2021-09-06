from Contour2D.CGrid2D import CGrids2D
from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import SparseGridIntegrator


def intfunc1(x: complex, y: complex) -> complex:
    s = 1
    m1sq = 0.2
    m2sq = 0.1
    return (1 - x) / (-s * x * (1 - x) * y - m2sq * (-1 + x + (1 - x) * y) + m1sq * (x + (1 - x) * y))


def intfunc2(x: complex, y: complex) -> complex:
    s = 1
    m2sq = 0.1
    return (1 - x) / (-s * x * (1 - x) * y - m2sq * (-1 + x + (1 - x) * y))


intgrand1 = Integrand2D(intfunc1, 0, 1, 0, 1)
intgrator = SparseGridIntegrator()
grid2d = CGrids2D(7, 7, 1, intgrator, intgrand1)
print(grid2d.Integrate())
print(grid2d.GatherInfo())

intgrand2 = Integrand2D(intfunc2, 0, 1, 0, 1)
intgrator2 = SparseGridIntegrator(epsilon=1.0e-3)
grid2dbigger = CGrids2D(103, 103, 1, intgrator2, intgrand2)
print(grid2dbigger.Integrate())
print(grid2dbigger.GatherInfo())
# grid2dbigger.Show()
