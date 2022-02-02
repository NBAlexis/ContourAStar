from Contour4D.Integrand4D import Integrand4D
from Contour4D.Integrator4D import SparseGridIntegrator4D


def intfunc(x: complex, y: complex, z: complex, w: complex) -> complex:
    return 1 / (x + y * y + z + w * w * w)


intgrand = Integrand4D(intfunc, -1, 1, 0, 1, -1, 1, -1, 1)
intgrator = SparseGridIntegrator4D()
print(intgrator.Integrate(intgrand, -1j, 1, 0, 1, 0, 1, 0, 1))
print(intgrand.GetDebugInfo())
# print(intgrator.Integrate(intgrand, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j))
