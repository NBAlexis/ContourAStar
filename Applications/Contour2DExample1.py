from Contour1D.Integrators import *

from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import SparseGridIntegrator
from Contour2D.ExtendPath2D import ExtendPath2D


def intfunc1(x: complex, y: complex) -> complex:
    return 1 / ((x + y - 0.25) ** 2)


intgrand1 = Integrand2D(intfunc1, 0, 1, 0, 1)
intgrand1.SetMathematicaExpress("""1 / ((x + y - 0.25) ^ 2)""")


def intfunc2(x: complex, y: complex) -> complex:
    return cmath.log(x + y - 0.25) / ((x + y - 0.25) ** 2)


intgrand2 = Integrand2D(intfunc2, 0, 1, 0, 1)
intgrand2.SetMathematicaExpress("""Log[x + y - 0.25] / ((x + y - 0.25) ^ 2)""")


def intfunc3(x: complex, y: complex) -> complex:
    return 1 / ((x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) ** 4)


intgrand3 = Integrand2D(intfunc3, 0, 1, 0, 1)
intgrand3.SetMathematicaExpress(
    """1 / ((x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) ^ 4)""")


def intfunc4(x: complex, y: complex) -> complex:
    return cmath.log(x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) / \
           ((x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) ** 4)


intgrand4 = Integrand2D(intfunc4, 0, 1, 0, 1)
intgrand4.SetMathematicaExpress("""Log[x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5] 
/ ((x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) ^ 4)""")


def intfunc5(x: complex, y: complex) -> complex:
    B = 5
    return (4 * (x - 1) * (x - 1) * x * y + x - 1) / (1 + B * (x - 1) * x * y)


intgrand5 = Integrand2D(intfunc5, 0, 1, 0, 1)
intgrand5.SetMathematicaExpress("""(4 * (x - 1) * (x - 1) * x * y + x - 1) / (1 + 5 * (x - 1) * x * y)""")


def intfunc6(x: complex, y: complex) -> complex:
    B = 100
    return (4 * (x - 1) * (x - 1) * x * y + x - 1) / (1 + B * (x - 1) * x * y)


intgrand6 = Integrand2D(intfunc6, 0, 1, 0, 1)
intgrand6.SetMathematicaExpress("""(4 * (x - 1) * (x - 1) * x * y + x - 1) / (1 + 100 * (x - 1) * x * y)""")


def intfunc7(x: complex, y: complex) -> complex:
    B = 10000
    return (4 * (x - 1) * (x - 1) * x * y + x - 1) / (1 + B * (x - 1) * x * y)


intgrand7 = Integrand2D(intfunc7, 0, 1, 0, 1)
intgrand7.SetMathematicaExpress("""(4 * (x - 1) * (x - 1) * x * y + x - 1) / (1 + 10000 * (x - 1) * x * y)""")


def intfunc8(x: complex, y: complex) -> complex:
    s = 1
    m1sq = 0.2
    m2sq = 0.1
    return (1 - x) / (-s * x * (1 - x) * y - m2sq * (-1 + x + (1 - x) * y) + m1sq * (x + (1 - x) * y))


intgrand8 = Integrand2D(intfunc8, 0, 1, 0, 1)
intgrand8.SetMathematicaExpress(
    """(1 - x) / (-1 * x * (1 - x) * y - 0.1 * (-1 + x + (1 - x) * y) + 0.2 * (x + (1 - x) * y))""")


intgrator = SparseGridIntegrator(logLevel=LogLevel.Warning)
gridsmall = ExtendPath2D(3, 3, 0, intgrator, logLevel=LogLevel.Warning)
gridbig = ExtendPath2D(11, 11, 1, intgrator, logLevel=LogLevel.General)
print(gridsmall.Integrate(intgrand1))
print(gridsmall.GatherInfo())

print(gridsmall.Integrate(intgrand2))
print(gridsmall.GatherInfo())

print(gridsmall.Integrate(intgrand3))
print(gridsmall.GatherInfo())

print(gridsmall.Integrate(intgrand4))
print(gridsmall.GatherInfo())

print(gridsmall.Integrate(intgrand5))
print(gridsmall.GatherInfo())

print(gridsmall.Integrate(intgrand6))
print(gridsmall.GatherInfo())

print(gridbig.Integrate(intgrand7))
print(gridbig.GatherInfo())

print(gridsmall.Integrate(intgrand8))
print(gridsmall.GatherInfo())

