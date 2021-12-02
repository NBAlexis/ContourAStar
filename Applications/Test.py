import cmath

import scipy.integrate

from Contour1D.CGrids import CGrids
from Contour1D.Integrand import *
from Contour1D.Integrators import *

# def testintegrand(x: complex) -> complex:
#     return 1 / (cmath.log(x) - 0.5)


# print(scipy.integrate.quad(testintegrand, 0, 1))
from Contour2D.CGrid2D import CGrids2D
from Contour2D.DisconnectedGrids import DisconnectedGrids2D
from Contour2D.Integrand2D import Integrand2D, XoverOneXMapping, XoverOneXSquareMapping, LogMapping
from Contour2D.Integrator2D import SparseGrid, SparseGridIntegrator

# from UsefulFunctions.SparseGridGenerator import TestSparseGridPoints, GetSparseGridNewPoints, GetWeightGaussPatterson, \
#     TestGridWithWeight

# TestGridWithWeight(2)

# def intfunc(x: complex, y: complex) -> complex:
#     return 1 / (cmath.log(x + y) + 1.0)
# intgrand = Integrand2D(intfunc, -1, 1, -1, 1)
from Contour3D.CGrid3D import CGrids3D
from Contour3D.CMarkedGrid3D import CMarkedGrids3D
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator3D
from UsefulFunctions.ClenshawCurtis import TestClenshawCurtisWeightList, ClenshawCurtis
from UsefulFunctions.ClenshawCurtisExp import ClenshawCurtisExp, TestClenshawCurtisExpWeightList
from UsefulFunctions.GaussianPattersonGenerator import GenerateGaussianPattersonOneOrder
from UsefulFunctions.QuadratureCache import SaveQuadrature
from UsefulFunctions.GaussianPatterson import TestGaussPattersonWeightList, GaussPatterson
from UsefulFunctions.NestedQuadrature import TestTrapezoidalWeightList, Trapezoidal
from UsefulFunctions.SparseGridGenerator import TestSparseGridPoints, TestGridWithWeight
from UsefulFunctions.SparseGridGenerator3D import TestSparseGridPoints3D, TestGridWithWeight3D


def intfuncB(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)


# """
# - 0.3 * y * x * x this term
def intfunc(x: complex, y: complex, z: complex) -> complex:
    return (1 - x) * (1 - x) * (1 - y) / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
                                           + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                                           + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                                           + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                                           - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                                           + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                                           + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                                           + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                                           + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2)


def intfuncC(x: complex, y: complex, z: complex) -> complex:
    return 1 / (abs(x + 1) ** 2 + abs(y + 1) ** 2 + abs(z + 1) ** 2 - 4.1)


def intfuncA(x: complex, y: complex, z: complex) -> complex:
    return 1 / (cmath.log((x + y + z)) + 0.5)


# """

intgrand = Integrand3D(intfuncC, -1, 1, -1, 1, -1, 1)
intgrand.SetMathematicaExpress("""1 / ((Abs[x + 1]^2 + Abs[y + 1]^2 + Abs[z + 1]^2 - 7.9))""")

intgrator = MathLinkIntegrator3D()
# print(intgrator.Integrate(intgrand, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j))
# print(intgrator.IntegrateYZ(intgrand, -1 - 1j, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j))
# print(intgrator.IntegrateZ(intgrand, -1 - 1j, 1 + 1j, -1 - 1j, 1 + 1j))

# grid3d = CMarkedGrids3D(3, 3, 0, intgrator, intgrand, requireEdge=False)
grid3d = CMarkedGrids3D(5, 5, 1, intgrator, intgrand, requireEdge=True)
print(grid3d.Integrate())
print(grid3d.GatherInfo())

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
