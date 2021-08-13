import cmath

import scipy.integrate

from Contour1D.CGrids import CGrids
from Contour1D.Integrand import *
from Contour1D.Integrators import *

# def testintegrand(x: complex) -> complex:
#     return 1 / (cmath.log(x) - 0.5)


# print(scipy.integrate.quad(testintegrand, 0, 1))
from Contour2D.CGrid2D import CGrids2D
from Contour2D.Integrand2D import Integrand2D, XoverOneXMapping, XoverOneXSquareMapping, LogMapping
from Contour2D.Integrator2D import SparseGrid, SparseGridIntegrator


# from UsefulFunctions.SparseGridGenerator import TestSparseGridPoints, GetSparseGridNewPoints, GetWeightGaussPatterson, \
#     TestGridWithWeight

# TestGridWithWeight(2)

# def intfunc(x: complex, y: complex) -> complex:
#     return 1 / (cmath.log(x + y) + 1.0)
# intgrand = Integrand2D(intfunc, -1, 1, -1, 1)
from UsefulFunctions.GaussianPatterson import TestGaussPattersonWeightList, GaussPatterson
from UsefulFunctions.NestedQuadrature import TestTrapezoidalWeightList, Trapezoidal
from UsefulFunctions.SparseGridGenerator import TestSparseGridPoints, TestGridWithWeight


def intfunc(y: complex, z: complex) -> complex:
    X = 4 + 1j
    return (y / (y * y + 0.25)) * cmath.sqrt(1 - z * z) * (-X + 2 * cmath.sqrt(X) * cmath.sqrt(y) * z) \
           / (X + y - 2 * cmath.sqrt(X) * cmath.sqrt(y) * z - 0.5j)

def intfunc2(x: complex, y: complex) -> complex:
    B = 10000
    return (4 * (x - 1) * (x - 1) * x * y + x - 1) \
           / (1 + B * (x - 1) * x * y)

# intgrand = Integrand2D(intfunc2, 0.3, 1, 0, 1)
intgrand = Integrand2D(intfunc2, 0, 1, 0, 1)
# intgrand = Integrand2D(intfunc, 0, cmath.inf, -1, 1, aInfMapping=LogMapping())
# intgrator = SparseGridIntegrator(nestedQuadrature=Trapezoidal(), maxOrder=11)
intgrator = SparseGridIntegrator(nestedQuadrature=GaussPatterson())

# print(intgrand.GetDebugInfo())
# print(intgrator.PartialIntegrateYEdge(intgrand, True, -1, 1))
# print(intgrator.PartialIntegrateYEdge(intgrand, False, -1, 1))
# print(intgrator.Integrate(intgrand, -1, -1 - 3j, -1, 1))
# print(intgrator.Integrate(intgrand, -1, -0.5, 0, -0.5j))
# print(intgrator.Integrate(intgrand, -1, -0.5, -0.5j, 1-0.5j))
# print(intgrator.Integrate(intgrand, -1, -0.5, 1-0.5j, 1))
# TestTrapezoidalWeightList()
# TestGaussPattersonWeightList()
grid2d = CGrids2D(31, 31, 5, intgrator, intgrand)
print(grid2d.Integrate())
print(grid2d.GatherInfo())
# grid2d.Show()

# trape = Trapezoidal()
# TestGridWithWeight(2, trape)
# TestGaussPattersonWeightList()