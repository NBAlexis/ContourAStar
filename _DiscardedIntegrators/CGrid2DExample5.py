from Contour1D.Integrators import *

from _DiscardedIntegrators.CGrid2D import CGrids2D
from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import SparseGridIntegrator


def intfunc1(x: complex, y: complex) -> complex:
    return 1 / ((x + y - 0.25) ** 2)


def intfunc2(x: complex, y: complex) -> complex:
    return cmath.log(x + y - 0.25) / ((x + y - 0.25) ** 2)


def intfunc3(x: complex, y: complex) -> complex:
    return 1 / ((x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) ** 4)


def intfunc4(x: complex, y: complex) -> complex:
    return cmath.log(x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) / (
                (x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) ** 4)


intgrand1 = Integrand2D(intfunc1, 0, 1, 0, 1)
intgrand2 = Integrand2D(intfunc2, 0, 1, 0, 1)
intgrand3 = Integrand2D(intfunc2, 0, 1, 0, 1)
intgrand4 = Integrand2D(intfunc2, 0, 1, 0, 1)
intgrator = SparseGridIntegrator()
grid2d = CGrids2D(3, 3, 0, intgrator, intgrand1)
print(grid2d.Integrate())
print(grid2d.GatherInfo())

grid2d.SetIntegrand(intgrand2)
grid2d.ResetGrid()
print(grid2d.Integrate())
print(grid2d.GatherInfo())

grid2d.SetIntegrand(intgrand3)
grid2d.ResetGrid()
print(grid2d.Integrate())
print(grid2d.GatherInfo())

grid2d.SetIntegrand(intgrand4)
grid2d.ResetGrid()
print(grid2d.Integrate())
print(grid2d.GatherInfo())