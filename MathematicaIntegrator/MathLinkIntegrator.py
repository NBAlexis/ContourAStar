import cmath

from Contour1D.CommonDefinitions import LogLevel
from Contour2D.Integrand2D import Integrand2D
from Contour2D.Integrator2D import Integrators2D
from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import Integrators3D
from MathematicaIntegrator.Constants import MathLink

_FAIL_WARNINGS_ = [
    "NIntegrate::slwcon",
    "NIntegrate::ncvb",
    "General::infy",
    "General::indet",
    "NIntegrate::inumri",
    "NIntegrate::inumr",
    "NIntegrate::eincr",
    "NIntegrate::oidiv",
    "NIntegrate::errprec"
]

_QUIET_WARNING_ = [
    "General::munfl",
    "NIntegrate::izero",
    "General::precw",

    "NIntegrate::slwcon",
    "NIntegrate::ncvb",
    "General::infy",
    "General::indet",
    "NIntegrate::inumri",
    "NIntegrate::inumr",
    "NIntegrate::eincr",
    "NIntegrate::oidiv",
    "NIntegrate::errprec"
]

_NINTEGRATE_CMD_ = 'Quiet[If[res = Check[NIntegrate[_integrand_, _range_ _op_], False, _warning1_]; NumberQ[res], {True, Re[res], Im[res]}, {False, 0, 0}], _warning2_]'


class MathLinkIntegrator2D(Integrators2D):

    def __init__(self, logLevel=LogLevel.Verbose):
        self.logLevel = logLevel
        self.mathlink = None
        self.lastFunc = None
        self.options = ""
        self.warningList = None
        self.quietList = None

    def __del__(self):
        print("========= MathLink Quit ========")
        self.Finish()

    def Integrate(self, func: Integrand2D,
                  fromX: complex, toX: complex,
                  fromY: complex, toY: complex) -> [bool, complex]:
        if not self.CheckFunc(func):
            return [False, cmath.nan]
        funcString = "g[x, y]"
        rangeString = "{" + "x, {}, {}".format(str(fromX).replace("j", " I"), str(toX).replace("j", " I")) + "}, " \
                      + "{" + "y, {}, {}".format(str(fromY).replace("j", " I"), str(toY).replace("j", " I")) + "}"
        return self.RealIntegrate(funcString, rangeString)

    def PartialIntegrateY(self, func: Integrand2D, x: complex, fromY: complex, toY: complex) -> [bool, complex]:
        if not self.CheckFunc(func):
            return [False, cmath.nan]
        funcString = "g[{}, y]".format(str(x).replace("j", " I"))
        rangeString = "{" + "z, {}, {}".format(str(fromY).replace("j", " I"), str(toY).replace("j", " I")) + "}"
        return self.RealIntegrate(funcString, rangeString)

    def CheckFunc(self, func: Integrand2D) -> bool:
        [hasFunc, funcExpression] = func.GetMathematicaExpress()
        if not hasFunc:
            return False
        if func != self.lastFunc:
            if self.mathlink is not None:
                self.mathlink.Quit()
                self.mathlink = None
        if self.mathlink is None:
            self.mathlink = MathLink()
            self.mathlink.Call(funcExpression)
            if self.logLevel >= LogLevel.Verbose:
                print("Put integrate function: " + funcExpression)
        self.lastFunc = func
        return True

    def SetOption(self, op: str):
        self.options = op

    def SetWarningList(self, warningList: list, quietList: list):
        self.warningList = warningList
        self.quietList = quietList

    def RealIntegrate(self, funcString: str, rangeString: str) -> [bool, complex]:
        failWarning = _FAIL_WARNINGS_ if self.warningList is None else self.warningList
        quietWarning = _QUIET_WARNING_ if self.quietList is None else self.quietList
        warningString1 = "{"
        for warningItem in failWarning:
            warningString1 = warningString1 + warningItem + ", "
        warningString1 = warningString1[0: len(warningString1) - 2] + "}"
        warningString2 = "{"
        for warningItem in quietWarning:
            warningString2 = warningString2 + warningItem + ", "
        warningString2 = warningString2[0: len(warningString2) - 2] + "}"
        cmdString = _NINTEGRATE_CMD_
        cmdString = cmdString.replace("_integrand_", funcString)
        cmdString = cmdString.replace("_range_", rangeString)
        cmdString = cmdString.replace("_warning1_", warningString1)
        cmdString = cmdString.replace("_warning2_", warningString2)
        cmdString = cmdString.replace(" _op_", self.options)
        if self.logLevel >= LogLevel.Verbose:
            print(cmdString)
        [bDone, resReal, resImage] = self.mathlink.Call(cmdString)
        if self.logLevel >= LogLevel.Verbose:
            print(bDone, " ", resReal, " ", resImage)
        if bDone:
            return [True, float(resReal) + float(resImage) * 1j]
        return [False, cmath.nan]

    def GetLeftEdgeX(self) -> float:
        return -1

    def Finish(self):
        if self.mathlink is not None:
            self.mathlink.Quit()
            self.mathlink = None
            self.lastFunc = None


class MathLinkIntegrator3D(Integrators3D):

    def __init__(self, logLevel=LogLevel.Verbose):
        self.logLevel = logLevel
        self.mathlink = None
        self.lastFunc = None
        self.options = ""
        self.warningList = None
        self.quietList = None

    def __del__(self):
        print("========= MathLink Quit ========")
        self.Finish()

    def Integrate(self, func: Integrand3D,
                  fromX: complex, toX: complex,
                  fromY: complex, toY: complex,
                  fromZ: complex, toZ: complex) -> [bool, complex]:
        if not self.CheckFunc(func):
            return [False, cmath.nan]
        funcString = "g[x, y, z]"
        rangeString = "{" + "x, {}, {}".format(str(fromX).replace("j", " I"), str(toX).replace("j", " I")) + "}, " \
                      + "{" + "y, {}, {}".format(str(fromY).replace("j", " I"), str(toY).replace("j", " I")) + "}, " \
                      + "{" + "z, {}, {}".format(str(fromZ).replace("j", " I"), str(toZ).replace("j", " I")) + "}"
        return self.RealIntegrate(funcString, rangeString)

    def IntegrateYZ(self, func: Integrand3D,
                    x: complex,
                    fromY: complex, toY: complex,
                    fromZ: complex, toZ: complex) -> [bool, complex]:
        if not self.CheckFunc(func):
            return [False, cmath.nan]
        funcString = "g[{}, y, z]".format(str(x).replace("j", " I"))
        rangeString = "{" + "y, {}, {}".format(str(fromY).replace("j", " I"), str(toY).replace("j", " I")) + "}, " \
                      + "{" + "z, {}, {}".format(str(fromZ).replace("j", " I"), str(toZ).replace("j", " I")) + "}"
        return self.RealIntegrate(funcString, rangeString)

    def IntegrateZ(self, func: Integrand3D,
                   x: complex, y: complex,
                   fromZ: complex, toZ: complex) -> [bool, complex]:
        if not self.CheckFunc(func):
            return [False, cmath.nan]
        funcString = "g[{}, {}, z]".format(str(x).replace("j", " I"), str(y).replace("j", " I"))
        rangeString = "{" + "z, {}, {}".format(str(fromZ).replace("j", " I"), str(toZ).replace("j", " I")) + "}"
        return self.RealIntegrate(funcString, rangeString)

    def CheckFunc(self, func: Integrand3D) -> bool:
        [hasFunc, funcExpression] = func.GetMathematicaExpress()
        if not hasFunc:
            return False
        if func != self.lastFunc:
            if self.mathlink is not None:
                self.mathlink.Quit()
                self.mathlink = None
        if self.mathlink is None:
            self.mathlink = MathLink()
            self.mathlink.Call(funcExpression)
            if self.logLevel >= LogLevel.Verbose:
                print("Put integrate function: " + funcExpression)
        self.lastFunc = func
        return True

    def SetOption(self, op: str):
        self.options = op

    def SetWarningList(self, warningList: list, quietList: list):
        self.warningList = warningList
        self.quietList = quietList

    def RealIntegrate(self, funcString: str, rangeString: str) -> [bool, complex]:
        failWarning = _FAIL_WARNINGS_ if self.warningList is None else self.warningList
        quietWarning = _QUIET_WARNING_ if self.quietList is None else self.quietList
        warningString1 = "{"
        for warningItem in failWarning:
            warningString1 = warningString1 + warningItem + ", "
        warningString1 = warningString1[0: len(warningString1) - 2] + "}"
        warningString2 = "{"
        for warningItem in quietWarning:
            warningString2 = warningString2 + warningItem + ", "
        warningString2 = warningString2[0: len(warningString2) - 2] + "}"
        cmdString = _NINTEGRATE_CMD_
        cmdString = cmdString.replace("_integrand_", funcString)
        cmdString = cmdString.replace("_range_", rangeString)
        cmdString = cmdString.replace("_warning1_", warningString1)
        cmdString = cmdString.replace("_warning2_", warningString2)
        cmdString = cmdString.replace(" _op_", self.options)
        if self.logLevel >= LogLevel.Verbose:
            print(cmdString)
        [bDone, resReal, resImage] = self.mathlink.Call(cmdString)
        if self.logLevel >= LogLevel.Verbose:
            print(bDone, " ", resReal, " ", resImage)
        if bDone:
            return [True, float(resReal) + float(resImage) * 1j]
        return [False, cmath.nan]

    def GetLeftEdgeX(self) -> float:
        return -1

    def Finish(self):
        if self.mathlink is not None:
            self.mathlink.Quit()
            self.mathlink = None
            self.lastFunc = None
