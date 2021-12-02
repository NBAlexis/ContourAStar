from Contour3D.CMarkedGrid3D import CMarkedGrids3D
from Contour3D.Integrand3D import Integrand3D
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator3D


def intfunc(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ** 2)


intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
intgrand.SetMathematicaExpress("""1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ^ 2)""")

intgrator = MathLinkIntegrator3D()
grid3d = CMarkedGrids3D(7, 7, 1, intgrator, intgrand, requireEdge=True)
print(grid3d.Integrate())
print(grid3d.GatherInfo())
intgrator.Finish()
