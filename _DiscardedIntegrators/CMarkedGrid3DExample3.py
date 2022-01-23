from _DiscardedIntegrators.CMarkedGrid3D import CMarkedGrids3D
from Contour3D.Integrand3D import Integrand3D
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator3D


def intfunc(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((-0.02 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)


intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
intgrand.SetMathematicaExpress("""(1 - x) * (1 - x) * (1 - y) / ((-0.02 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ^ 2)""")

intgrator = MathLinkIntegrator3D()
grid3d = CMarkedGrids3D(7, 7, 1, intgrator, intgrand, requireEdge=True)
print(grid3d.Integrate())
print(grid3d.GatherInfo())
intgrator.Finish()
