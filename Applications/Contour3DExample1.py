import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour3D.ExtendPath3D import ExtendPath3D
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D


def intfunc1(x: complex, y: complex, z: complex) -> complex:
    return 1 / (cmath.log((x + y + z)) + 0.5)


def intfunc1d(x: complex, y: complex, z: complex) -> complex:
    return cmath.log((x + y + z)) + 0.5


intgrand1 = Integrand3D(intfunc1, -1, 1, -1, 1, -1, 1)
intgrand1.SetMathematicaExpress("""1 / (Log[x + y + z] + 0.5)""")
intgrand1.SetDenominator(intfunc1d)


def intfunc2(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((x + y + z - 0.25) ** 2)


def intfunc2d(x: complex, y: complex, z: complex) -> complex:
    return (x + y + z - 0.25) ** 2


intgrand2 = Integrand3D(intfunc2, 0, 1, 0, 1, 0, 1)
intgrand2.SetMathematicaExpress("""1 / ((x + y + z - 0.25) ^ 2)""")
intgrand2.SetDenominator(intfunc2d)


def intfunc3(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((x + y + z + x * x - 0.25) ** 2)


def intfunc3d(x: complex, y: complex, z: complex) -> complex:
    return (x + y + z + x * x - 0.25) ** 2


intgrand3 = Integrand3D(intfunc3, 0, 1, 0, 1, 0, 1)
intgrand3.SetMathematicaExpress("""1 / ((x + y + z + x * x - 0.25) ^ 2)""")
intgrand3.SetDenominator(intfunc3d)


def intfunc4(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ** 2)


def intfunc4d(x: complex, y: complex, z: complex) -> complex:
    return (0.1 * x + 0.2 * y + 0.3 * z
            + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
            + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
            + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ** 2


intgrand4 = Integrand3D(intfunc4, 0, 1, 0, 1, 0, 1)
intgrand4.SetMathematicaExpress("""1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ^ 2)""")
intgrand4.SetDenominator(intfunc4d)


def intfunc5(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2)


def intfunc5d(x: complex, y: complex, z: complex) -> complex:
    return (0.1 * x + 0.2 * y + 0.3 * z
            + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
            + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
            + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
            - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
            + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
            + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
            + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
            + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2


intgrand5 = Integrand3D(intfunc5, 0, 1, 0, 1, 0, 1)
intgrand5.SetMathematicaExpress("""1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2)""")
intgrand5.SetDenominator(intfunc5d)


def intfunc6(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2)


def intfunc6d(x: complex, y: complex, z: complex) -> complex:
    return (0.1 * x + 0.2 * y + 0.3 * z
            + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
            + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
            + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
            - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
            + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
            + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
            + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
            - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2


intgrand6 = Integrand3D(intfunc6, 0, 1, 0, 1, 0, 1)
intgrand6.SetMathematicaExpress("""1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2)""")
intgrand6.SetDenominator(intfunc6d)


def intfunc7(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((-0.02 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)


def intfunc7d(x: complex, y: complex, z: complex) -> complex:
    return (-0.02 - (1 - x) * (1 - y) * x * z
            + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2


intgrand7 = Integrand3D(intfunc7, 0, 1, 0, 1, 0, 1)
intgrand7.SetMathematicaExpress("""(1 - x) * (1 - x) * (1 - y) / ((-0.02 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ^ 2)""")
intgrand7.SetDenominator(intfunc7d)


def intfunc8(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)


def intfunc8d(x: complex, y: complex, z: complex) -> complex:
    return (0.2 - (1 - x) * (1 - y) * x * z
            + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2


intgrand8 = Integrand3D(intfunc8, 0, 1, 0, 1, 0, 1)
intgrand8.SetMathematicaExpress("""(1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ^ 2)""")
intgrand8.SetDenominator(intfunc8d)


def intfunc9(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((0.2 - 0.142857 * x * (1 - x) * y
                                           - x * (1 - x) * (1 - y) * z
                                           - 0.142857 * (1 - x) * (1 - x) * y * (1 - y) * z
                                           + 0.4 * (1 - x) * (1 - x) * y * (1 - y) * (1 - z)) ** 2)


def intfunc9d(x: complex, y: complex, z: complex) -> complex:
    return (0.2 - 0.142857 * x * (1 - x) * y
            - x * (1 - x) * (1 - y) * z
            - 0.142857 * (1 - x) * (1 - x) * y * (1 - y) * z
            + 0.4 * (1 - x) * (1 - x) * y * (1 - y) * (1 - z)) ** 2


intgrand9 = Integrand3D(intfunc9, 0, 1, 0, 1, 0, 1)

intgrand9.SetMathematicaExpress("""(1-x)*(1-x)*(1-y)/((0.2 - 0.142857*x*(1-x)*y
                               - x*(1-x)*(1-y)*z
                               - 0.142857*(1-x)*(1-x)*y*(1-y)*z 
                               + 0.4*(1-x)*(1-x)*y*(1-y)*(1-z))^2)""")
intgrand9.SetDenominator(intfunc9d)

intgrator = SparseGridIntegrator3D(logLevel=LogLevel.Warning)
grid3dsmall = ExtendPath3D(3, 3, 0, intgrator, logLevel=LogLevel.Warning)
grid3d = ExtendPath3D(5, 5, 0, intgrator, logLevel=LogLevel.Warning)
grid3dmid = ExtendPath3D(13, 13, 1, intgrator, logLevel=LogLevel.Warning)

# This one is special, 5x5 grid will not work!!(why??)
print(grid3dsmall.Integrate(intgrand1))
print(grid3dsmall.GatherInfo())

print(grid3d.Integrate(intgrand2))
print(grid3d.GatherInfo())

print(grid3d.Integrate(intgrand3))
print(grid3d.GatherInfo())

print(grid3d.Integrate(intgrand4))
print(grid3d.GatherInfo())

print(grid3d.Integrate(intgrand5))
print(grid3d.GatherInfo())

print(grid3d.Integrate(intgrand6))
print(grid3d.GatherInfo())

print(grid3dmid.Integrate(intgrand7))
print(grid3dmid.GatherInfo())

print(grid3d.Integrate(intgrand8))
print(grid3d.GatherInfo())

print(grid3d.Integrate(intgrand9))
print(grid3d.GatherInfo())
