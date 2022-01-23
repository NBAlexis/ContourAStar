from datetime import datetime

import numpy as np

from Contour1D.CommonDefinitions import LogLevel
from Contour3D.Integrand3D import Integrand3D
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator3D
from MonteCarlo.MCIntegrator3D import MCIntegrator3D


def intfunc(x: complex, y: complex, z: complex) -> complex:
    return 1.0 / ((3.0001 - x ** 2 - y ** 2 - z ** 2) ** 3)


vpf = np.vectorize(intfunc)

intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
intgrand.MakeSureZeroOne()
intgrand.SetMathematicaExpress("""1 / ((30001/10000 - x^2 - y^2 - z^2) ^ 3)""")
intgrand.SetVPF(vpf)

print(datetime.now())
intgrator1 = MCIntegrator3D(pointCount=10000)
print(datetime.now())
intgratorM = MathLinkIntegrator3D(logLevel=LogLevel.Verbose)

print(datetime.now())
print(intgrator1.Integrate(intgrand, 0, 1, 0, 1, 0, 1))
print(datetime.now())
intgratorM.Integrate(intgrand, 0, 1, 0, 1, 0, 1)
intgratorM.Finish()
