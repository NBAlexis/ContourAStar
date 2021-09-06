import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour3D.CGrid3D import CGrids3D
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D


def intfunc(x: complex, y: complex, z: complex) -> complex:
    return 1 / (cmath.log((x + y + z)) + 0.5)


intgrand = Integrand3D(intfunc, -1, 1, -1, 1, -1, 1)
intgrator = SparseGridIntegrator3D(logLevel=LogLevel.Verbose)
grid3d = CGrids3D(3, 3, 0, intgrator, intgrand)
print(grid3d.Integrate())
print(grid3d.GatherInfo())

"""
(-3*ExpIntegralEi[3/2] + 3*ExpIntegralEi[3/2 + (3*I)*Pi] + 3*E*(-ExpIntegralEi[1/2] + ExpIntegralEi[1/2 + I*Pi] + 3*ExpIntegralEi[1/2 + Log[3]] - 3*ExpIntegralEi[1/2 + I*Pi + Log[3]]) + 
  6*Sqrt[E]*(ExpIntegralEi[1] + ExpIntegralEi[1 + (2*I)*Pi] - ExpIntegralEi[1 + Log[9]] - ExpIntegralEi[1 + (2*I)*Pi + Log[9]]) + ExpIntegralEi[3/2 + Log[27]] - 
  ExpIntegralEi[3/2 + (3*I)*Pi + Log[27]])/(2*E^(3/2))
"""