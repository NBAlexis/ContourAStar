from _DiscardedIntegrators.CGrid2D import CGrids2D
from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import SparseGridIntegrator


def intfunc(x: complex, y: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / (
                (0.2 + 0.4 * (1 - x) * (1 - x) * (1 - y) * y) * (0.2 - x * (1 - x - y + x * y)))


intgrand = Integrand2D(intfunc, 0, 1, 0, 1)
intgrator = SparseGridIntegrator()
grid2d = CGrids2D(3, 3, 0, intgrator, intgrand)
print(grid2d.Integrate())
print(grid2d.GatherInfo())

