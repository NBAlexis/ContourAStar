import cmath

from Contour2D.CGrid2D import CGrids2D
from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import SparseGridIntegrator


def intfunc2(x: complex, y: complex) -> complex:
    return 1 / (cmath.log(x + y) + 0.5)


intgrand = Integrand2D(intfunc2, -1, 1, -1, 1)
intgrator = SparseGridIntegrator()
grid2d = CGrids2D(3, 3, 0, intgrator, intgrand)
print(grid2d.Integrate())
print(grid2d.GatherInfo())
