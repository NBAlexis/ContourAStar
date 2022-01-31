import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour3D.ExtendPath3D import ExtendPath3D
from Contour3D.Integrand3D import Integrand3D
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator3D

"""
These integrations cannot be solved using Sparse grid
"""


def intfunc1(x: complex, y: complex, z: complex) -> complex:
    return 1 / (cmath.log((x + y + z)) + 0.5)


intgrand1 = Integrand3D(intfunc1, -1, 1, -1, 1, -1, 1)
intgrand1.SetMathematicaExpress("""1 / (Log[x + y + z] + 0.5)""")


def intfunc2(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ** 2)


def intfunc2d(x: complex, y: complex, z: complex) -> complex:
    return (0.2 - (1 - x) * (1 - y) * x * z
            - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ** 2


intgrand2 = Integrand3D(intfunc2, 0, 1, 0, 1, 0, 1)
intgrand2.SetMathematicaExpress("""(1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ^ 2)""")
intgrand2.SetDenominator(intfunc2d)


def intfunc3(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ** 2)


def intfunc3d(x: complex, y: complex, z: complex) -> complex:
    return (0.2 - (1 - x) * (1 - y) * x * z
            + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ** 2


intgrand3 = Integrand3D(intfunc3, 0, 1, 0, 1, 0, 1)
intgrand3.SetMathematicaExpress("""(1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ^ 2)""")
intgrand3.SetDenominator(intfunc3d)

intgrator = MathLinkIntegrator3D(logLevel=LogLevel.Warning)
grid3dsmall = ExtendPath3D(5, 5, 0, intgrator, logLevel=LogLevel.General)
grid3d = ExtendPath3D(11, 11, 0, intgrator, logLevel=LogLevel.General)

print(grid3dsmall.Integrate(intgrand1))
res1 = grid3dsmall.GatherInfo()

print(grid3d.Integrate(intgrand2))
res2 = grid3d.GatherInfo()

print(grid3d.Integrate(intgrand3))
res3 = grid3d.GatherInfo()

print("===res===")
print(res1)
print(res2)
print(res3)

intgrator.Finish()
