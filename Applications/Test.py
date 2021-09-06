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
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D
from UsefulFunctions.GaussianPatterson import TestGaussPattersonWeightList, GaussPatterson
from UsefulFunctions.NestedQuadrature import TestTrapezoidalWeightList, Trapezoidal
from UsefulFunctions.SparseGridGenerator import TestSparseGridPoints, TestGridWithWeight
from UsefulFunctions.SparseGridGenerator3D import TestSparseGridPoints3D, TestGridWithWeight3D

"""
# def intfunc(x: complex, y: complex, z: complex) -> complex:
#    return (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
#                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)
def intfunc(x: complex, y: complex, z: complex) -> complex:
     return 1 / ((x + y + z + x * x - 0.25) ** 2)


# TestGridWithWeight3D(3)
intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
intgrator = SparseGridIntegrator3D(epsilon=0.001)
grid3d = CGrids3D(33, 33, 11, intgrator, intgrand)
print(grid3d.Integrate())
print(grid3d.GatherInfo())
# grid3d.Show()
"""


# """


def intfunc(x: complex, y: complex) -> complex:
    return 1 / ((x + x * x + x * x * x + x * y + x * x * y + x * x * x * y + y * y * y + y * y + y - 0.5) ** 4)


intgrand = Integrand2D(intfunc, 0, 1, 0, 1)
intgrator = SparseGridIntegrator()
grid2d = DisconnectedGrids2D(7, 1, intgrator, intgrand, maxStep=10000000)
print(grid2d.Integrate())
# print(grid2d.GatherInfo())
# """
