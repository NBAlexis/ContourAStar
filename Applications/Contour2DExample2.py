from Contour2D.CGrid2D import CGrids2D
from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import SparseGridIntegrator


def intfunc1(x: complex, y: complex) -> complex:
    B = 5
    return (4 * (x - 1) * (x - 1) * x * y + x - 1) \
           / (1 + B * (x - 1) * x * y)


def intfunc2(x: complex, y: complex) -> complex:
    B = 100
    return (4 * (x - 1) * (x - 1) * x * y + x - 1) \
           / (1 + B * (x - 1) * x * y)


def intfunc3(x: complex, y: complex) -> complex:
    B = 10000
    return (4 * (x - 1) * (x - 1) * x * y + x - 1) \
           / (1 + B * (x - 1) * x * y)

intgrand = Integrand2D(intfunc1, 0, 1, 0, 1)
intgrator = SparseGridIntegrator()
grid2d = CGrids2D(7, 7, 0, intgrator, intgrand)
print(grid2d.Integrate())
print(grid2d.GatherInfo())

intgrand2 = Integrand2D(intfunc2, 0, 1, 0, 1)
grid2d.SetIntegrand(intgrand2)
grid2d.ResetGrid()
print(grid2d.Integrate())
print(grid2d.GatherInfo())

intgrand3 = Integrand2D(intfunc3, 0, 1, 0, 1)
grid2dbigger = CGrids2D(23, 23, 1, intgrator, intgrand3)
print(grid2dbigger.Integrate())
print(grid2dbigger.GatherInfo())

