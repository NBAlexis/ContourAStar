"""
In case there are poles at end-point
"""
from Contour1D.CommonDefinitions import LogLevel
from Contour2D.ExtendPath2D import ExtendPath2D
from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import SparseGridIntegrator
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator2D


def intfunc(x: complex, y: complex) -> complex:
    s = 1
    m2sq = 0.1
    return (1 - x) / (-s * x * (1 - x) * y - m2sq * (-1 + x + (1 - x) * y))


intgrand = Integrand2D(intfunc, 0, 1, 0, 1)
intgrand.SetMathematicaExpress("""(1 - x) / (-1 * x * (1 - x) * y - 0.1 * (-1 + x + (1 - x) * y))""")

intgrator = SparseGridIntegrator(logLevel=LogLevel.Warning, epsilon=2.0e-4)

gridhuge = ExtendPath2D(15, 7, 5, intgrator, logLevel=LogLevel.General)
print(gridhuge.Integrate(intgrand))
print(gridhuge.GatherInfo())

"""
If there was Mathematica installed
in the folder:

Mathematica/12.0/SystemFiles/Links/WolframClientForPython

use "pip install ."

Then, MathLink integration should work
"""


intgratorm = MathLinkIntegrator2D()

gridsmall = ExtendPath2D(3, 3, 0, intgratorm, logLevel=LogLevel.General)
print(gridsmall.Integrate(intgrand))
print(gridsmall.GatherInfo())

intgratorm.Finish()