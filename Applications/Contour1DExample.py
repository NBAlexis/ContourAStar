from Contour1D.CGrids import CGrids
from Contour1D.Integrand import *
from Contour1D.Integrators import *


def testIntegrand1(x: complex) -> complex:
    return 1 / (cmath.log(x - 0.5) + 0.5)


def testIntegrand2(x: complex) -> complex:
    return cmath.sqrt(1 - x * x) * (x - 0.25j) / (x - 0.5)


def testIntegrand3(x: complex) -> complex:
    return 1 / (cmath.exp(3 * x) - 5)


def testIntegrand4(x: complex) -> complex:
    return 1 / (cmath.log(x) + 0.5)


testIntegrator = Simpson(logLevel=LogLevel.Verbose)

"""
NIntegrate[1/(Log[a - 0.5] + 0.5), {a, 0, 1}] = -0.802623-0.136651 I
"""
integrand1 = Integrand(testIntegrand1, 0, 1, IntegrandType.ZeroOne)
grids1 = CGrids(3, 3, 0, testIntegrator, integrand1)
[res, resV] = grids1.Integrate()
print("{}: {}".format(res, resV))
grids1.ShowIntegralPath()

"""
NIntegrate[Sqrt[1 - z^2] ((z - I/4)/(z - 1/2)), {z, -1, -1 + I}] + 
 NIntegrate[Sqrt[1 - z^2] ((z - I/4)/(z - 1/2)), {z, -1 + I, 1}]
= 0.105223 -0.96765 I

NIntegrate[Sqrt[1 - z^2] ((z - I/4)/(z - 1/2)), {z, -1, -1 - I}] + 
 NIntegrate[Sqrt[1 - z^2] ((z - I/4)/(z - 1/2)), {z, -1 - I, 1}]
= 1.46557 +1.75305 I

see: ComplexPlot[Sqrt[1 - z^2] ((z - I/4)/(z - 1/2)), {z, -2 - 2 I, 2 + 2 I}]
"""
integrand2 = Integrand(testIntegrand2, -1, 1)
grids2 = CGrids(3, 3, 0, testIntegrator, integrand2)
[res, resV] = grids2.Integrate()
print("{}: {}".format(res, resV))
grids2.ShowIntegralPath()

"""
NIntegrate[1/(Exp[3 x] - 5), {x, 0.5, 0.5 + I}] + 
 NIntegrate[1/(Exp[3 x] - 5), {x, 0.5 + I, Infinity}]
= 0.143812 -0.20944 I
"""
integrand3 = Integrand(testIntegrand3, 0.5, cmath.inf)
grids3 = CGrids(3, 3, 0, testIntegrator, integrand3)
[res, resV] = grids3.Integrate()
print("{}: {}".format(res, resV))
grids3.ShowIntegralPath()


"""
NIntegrate[1/(Log[a] + 0.5), {a, 0, I}] + 
 NIntegrate[1/(Log[a] + 0.5), {a, I, 1 + I}] + 
 NIntegrate[1/(Log[a] + 0.5), {a, 1 + I, 1}]
= 0.275498 -1.90547 I
"""
integrand4 = Integrand(testIntegrand4, 0, 1, IntegrandType.ZeroOne)
grids4 = CGrids(3, 3, 0, testIntegrator, integrand4)
[res, resV] = grids4.Integrate()
print("{}: {}".format(res, resV))
grids4.ShowIntegralPath()
