from Contour3D.CGrid3D import CGrids3D
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D


def intfunc(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((x + y + z - 0.25) ** 2)


intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
intgrator = SparseGridIntegrator3D(epsilon=0.001)
grid3d = CGrids3D(7, 7, 1, intgrator, intgrand)
print(grid3d.Integrate())
print(grid3d.GatherInfo())

