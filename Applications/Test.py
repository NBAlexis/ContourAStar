# def testintegrand(x: complex) -> complex:
#     return 1 / (cmath.log(x) - 0.5)


# print(scipy.integrate.quad(testintegrand, 0, 1))

# from SparseGridIntegrators.SparseGridGenerator import TestSparseGridPoints, GetSparseGridNewPoints, GetWeightGaussPatterson, \
#     TestGridWithWeight

# TestGridWithWeight(2)

# def intfunc(x: complex, y: complex) -> complex:
#     return 1 / (cmath.log(x + y) + 1.0)
# intgrand = Integrand2D(intfunc, -1, 1, -1, 1)
import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour1D.IntegrandV2 import IntegrandV2
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D
from Contour3D.ExtendPath3D import ExtendPath3D
from Contour4D.ExtendPath4D import ExtendPath4D
from Contour4D.Integrand4D import Integrand4D
from Contour4D.Integrator4D import SparseGridIntegrator4D
from Contour4D.Integrator4DV2 import SparseGridIntegrator4DV2
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator3D


def intfuncA(x: complex, y: complex, z: complex, w: complex) -> complex:
    return 1 / (cmath.log(x + y + z + w) + 1)


def intfunc(x: complex, y: complex, z: complex, w: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((0.2 - 0.142857 * x * (1 - x) * y
                                           - x * (1 - x) * (1 - y) * z
                                           - 0.142857 * (1 - x) * (1 - x) * y * (1 - y) * z
                                           + 0.4 * (1 - x) * (1 - x) * y * (1 - y) * (1 - z)
                                           + w * w * x
                                           - w * w * y
                                           + z * z * w
                                           - w) ** 2)


def intfuncd(x: complex, y: complex, z: complex, w: complex) -> complex:
    return 0.2 - 0.142857 * x * (1 - x) * y \
           - x * (1 - x) * (1 - y) * z \
           - 0.142857 * (1 - x) * (1 - x) * y * (1 - y) * z \
           + 0.4 * (1 - x) * (1 - x) * y * (1 - y) * (1 - z) \
           + w * w * x \
           - w * w * y \
           + z * z * w \
           - w


# """

intgrand = IntegrandV2(intfunc, None, """(1 - x) * (1 - x) * (1 - y) / ((0.2 - 0.142857 * x * (1 - x) * y
                                           - x * (1 - x) * (1 - y) * z
                                           - 0.142857 * (1 - x) * (1 - x) * y * (1 - y) * z
                                           + 0.4 * (1 - x) * (1 - x) * y * (1 - y) * (1 - z)
                                           + w * w * x
                                           - w * w * y
                                           + z * z * w
                                           - w) ^ 2)""", "f[(w+1)/2,(w+1)/2,(w+1)/2,(w+1)/2]/16")

# intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
# intgrand.SetMathematicaExpress("""1 / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
#                  + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
#                  + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
#                  + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
#                  - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
#                  + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
#                  + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
#                  + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
#                  - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2)""")
# intgrator = MathLinkIntegrator3D()
intgrator = SparseGridIntegrator4DV2(epsilon=1.0e-2, logLevel=LogLevel.Warning)
# print(intgrator.Integrate(intgrand, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j))
# print(intgrator.IntegrateYZ(intgrand, -1 - 1j, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j))
# print(intgrator.IntegrateZ(intgrand, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j))

# grid3d = ExtendPath3D(11, 11, 0, intgrator, intgrand)
grid4d = ExtendPath4D(7, 7, 0, intgrator, intgrand)
# grid3d = CMarkedGrids3D(5, 5, 1, intgrator, intgrand, requireEdge=True)
print(grid4d.Integrate())
print(grid4d.GatherInfo())

intgrator.Finish()

"""
# TestGridWithWeight3D(3)
intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
intgrator = SparseGridIntegrator3D(epsilon=0.01, epsilon2d=0.1, epsilon1d=0.1, nestedQuadrature=ClenshawCurtis(), maxOrder=30)
grid3d = CMarkedGrids3D(23, 23, 1, intgrator, intgrand, requireEdge=False)
print(grid3d.Integrate())
print(grid3d.GatherInfo())
# grid3d.Show()
"""

# TestClenshawCurtisWeightList()
# TestClenshawCurtisWeightList()
# TestSparseGridPoints(6, ClenshawCurtis())
# TestClenshawCurtisWeightList()
# TestGridWithWeight(1, nestedQuadrature=ClenshawCurtis())

# TestClenshawCurtisExpWeightList()
# """
# def func(x, y):
#    return cmath.exp(x-y)

# intgrand = Integrand2D(func, -1, 1, -1, 1)
# intgrator = SparseGridIntegrator(nestedQuadrature=ClenshawCurtisExp(), maxOrder=15)
# print(intgrator.Integrate(intgrand, -1, 1, -1, 1))
# """

# SaveQuadrature("../_Data/ClenshawCurtisExp/", ClenshawCurtisExp(), 16)
# [x, w] = GenerateGaussianPattersonOneOrder(3)
# print(x)
# print(w)
