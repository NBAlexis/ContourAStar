import numpy as np

from Contour1D.CommonDefinitions import LogLevel
from Contour1D.IntegrandV2 import IntegrandV2
from Contour4D.ExtendPath4D import ExtendPath4D
from Contour4D.Integrator4DV2 import SparseGridIntegrator4DV2


def intfunc1(x: complex, y: complex, z: complex, w: complex) -> complex:
    return 1 / (np.log(x + y + z + w) + 1)


def intfunc1d(x: complex, y: complex, z: complex, w: complex) -> complex:
    return np.log(x + y + z + w) + 1


intgrand1 = IntegrandV2(intfunc1, intfunc1d, "1 / (Log[x+y+z+w]+1)", "f[x,y,z,w]")


intgrator = SparseGridIntegrator4DV2(logLevel=LogLevel.Warning)

grid4d = ExtendPath4D(3, 3, 0, intgrator, intgrand1)
print(grid4d.Integrate())
print(grid4d.GatherInfo())

intgrator.Finish()
