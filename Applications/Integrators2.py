from datetime import datetime

from Contour1D.CommonDefinitions import LogLevel
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D
from MathematicaIntegrator.MathLinkIntegrator import MathLinkIntegrator3D
from SparseGridIntegrators.ClenshawCurtis import ClenshawCurtis
from SparseGridIntegrators.ClenshawCurtisExp import ClenshawCurtisExp
from SparseGridIntegrators.NestedQuadrature import Trapezoidal


def intfunc(x: complex, y: complex, z: complex) -> complex:
    return 1.0 / ((3.0001 - x ** 2 - y ** 2 - z ** 2) ** 3)


intgrand = Integrand3D(intfunc, 0, 1, 0, 1, 0, 1)
intgrand.MakeSureZeroOne()
intgrand.SetMathematicaExpress("""1 / ((30001/10000 - x^2 - y^2 - z^2) ^ 3)""")

print(datetime.now())
intgrator1 = SparseGridIntegrator3D(logLevel=LogLevel.Verbose)
print(datetime.now())
intgrator2 = SparseGridIntegrator3D(nestedQuadrature=ClenshawCurtis())
print(datetime.now())
intgrator3 = SparseGridIntegrator3D(nestedQuadrature=ClenshawCurtisExp(), maxOrder=12)
print(datetime.now())
intgrator4 = SparseGridIntegrator3D(nestedQuadrature=Trapezoidal(), maxOrder=7)
print(datetime.now())
intgratorM = MathLinkIntegrator3D(logLevel=LogLevel.Verbose)
# intgratorM.SetWarningList([
#     "NIntegrate::ncvb",
#     "General::infy",
#     "General::indet",
#     "NIntegrate::inumri",
#     "NIntegrate::inumr",
#     "NIntegrate::eincr",
#     "NIntegrate::oidiv",
#     "NIntegrate::errprec"
# ], None)
# intgrator2.SetOption(", WorkingPrecision -> 10")

print(intgrator1.Integrate(intgrand, 0, 1, 0, 1, 0, 1))
print(intgrator2.Integrate(intgrand, 0, 1, 0, 1, 0, 1))
print(intgrator3.Integrate(intgrand, 0, 1, 0, 1, 0, 1))
print(intgrator4.Integrate(intgrand, 0, 1, 0, 1, 0, 1))
intgratorM.Integrate(intgrand, 0, 1, 0, 1, 0, 1)
intgratorM.Finish()
