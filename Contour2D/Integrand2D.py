import cmath

from Contour1D.Integrand import IntegrandType


class AInfiniteMapping:
    def Factor(self, bPM: bool, a: complex, x: complex) -> complex:
        pass

    def Arg(self, bPM: bool, a: complex, x: complex) -> complex:
        pass

    def FactorString(self, bPM: bool, a: complex, variableName: str) -> str:
        pass

    def ArgString(self, bPM: bool, a: complex, variableName: str) -> str:
        pass


class LogMapping(AInfiniteMapping):
    def Factor(self, bPM: bool, a: complex, x: complex) -> complex:
        return 1 / (x + 1)

    def Arg(self, bPM: bool, a: complex, x: complex) -> complex:
        return (a - cmath.log(0.5 * (x + 1))) if bPM else (a + cmath.log(0.5 * (x + 1)))

    def FactorString(self, bPM: bool, a: complex, variableName: str) -> str:
        return "(1 / ({} + 1)) * ".format(variableName)

    def ArgString(self, bPM: bool, a: complex, variableName: str) -> str:
        return "{} {} Log[({} + 1) / 2]".format(a, "-" if bPM else "+", variableName)


class XoverOneXMapping(AInfiniteMapping):
    def Factor(self, bPM: bool, a: complex, x: complex) -> complex:
        return 2 / ((x - 1) * (x - 1))

    def Arg(self, bPM: bool, a: complex, x: complex) -> complex:
        return (a + (1 + x) / (1 - x)) if bPM else (a - (1 + x) / (1 - x))

    def FactorString(self, bPM: bool, a: complex, variableName: str) -> str:
        return "2/(({}-1)^2) * ".format(variableName)

    def ArgString(self, bPM: bool, a: complex, variableName: str) -> str:
        return "{} {} (1+{})/(1-{})".format(a, "+" if bPM else "-", variableName, variableName)


class XoverOneXSquareMapping(AInfiniteMapping):
    def Factor(self, bPM: bool, a: complex, x: complex) -> complex:
        return 1 / ((x - 1) * (x - 1) * cmath.sqrt((1 + x) / (1 - x)))

    def Arg(self, bPM: bool, a: complex, x: complex) -> complex:
        return (a + cmath.sqrt((1 + x) / (1 - x))) if bPM else (a - cmath.sqrt((1 + x) / (1 - x)))

    def FactorString(self, bPM: bool, a: complex, variableName: str) -> str:
        return "1/(({}-1)^2 * Sqrt[(1+{})/(1-{})]) * ".format(variableName, variableName, variableName)

    def ArgString(self, bPM: bool, a: complex, variableName: str) -> str:
        return "{} {} Sqrt[(1+{})/(1-{})]".format(a, "+" if bPM else "-", variableName, variableName)


class Integrand2D:

    def __init__(self, func, leftX: complex, rightX: complex, leftY: complex, rightY: complex,
                 caseX: IntegrandType = IntegrandType.AB, caseY: IntegrandType = IntegrandType.AB,
                 aInfMapping: AInfiniteMapping = LogMapping()):
        self.func = func
        # ================= X ====================
        self.caseX: IntegrandType = caseX
        if cmath.isinf(leftX) and cmath.isinf(rightX):
            self.caseX = IntegrandType.InfiniteInfinite
        elif cmath.isinf(leftX):
            self.caseX = IntegrandType.InfiniteA
        elif cmath.isinf(rightX):
            self.caseX = IntegrandType.AInfinite
        self.leftX: complex = leftX
        self.rightX: complex = rightX
        self.sepX: complex = 0
        if caseX == IntegrandType.AB:
            self.sepX = (rightX - leftX) * 0.5
        # ================= Y ====================
        self.caseY: IntegrandType = caseY
        if cmath.isinf(leftY) and cmath.isinf(rightY):
            self.caseY = IntegrandType.InfiniteInfinite
        elif cmath.isinf(leftY):
            self.caseY = IntegrandType.InfiniteA
        elif cmath.isinf(rightY):
            self.caseY = IntegrandType.AInfinite
        self.leftY: complex = leftY
        self.rightY: complex = rightY
        self.sepY: complex = 0
        if caseY == IntegrandType.AB:
            self.sepY = (rightY - leftY) * 0.5
        self.aInfMapping = aInfMapping
        self.hasMathematicaExpression = False
        self.FExpression = ""

    def SetMathematicaExpress(self, expression: str) -> str:
        self.hasMathematicaExpression = True
        self.FExpression = "f[x_, y_]:=" + expression + ";\n" + self.Dress()
        return self.FExpression

    def GetMathematicaExpress(self) -> [bool, str]:
        if self.hasMathematicaExpression:
            return [True, self.FExpression]
        return [False, ""]

    def Dress(self) -> str:
        argX = ""
        argY = ""
        factor = ""
        if self.caseX == IntegrandType.AB:
            factor = "{} * ".format(self.sepX)
            argX = "{} * (x + 1) + {}".format(self.sepX, self.leftX)
        elif self.caseX == IntegrandType.AInfinite:
            factor = self.aInfMapping.FactorString(True, self.leftX, "x")
            argX = self.aInfMapping.ArgString(True, self.leftX, "x")
        elif self.caseX == IntegrandType.InfiniteA:
            factor = self.aInfMapping.FactorString(False, self.rightX, "x")
            argX = self.aInfMapping.ArgString(False, self.rightX, "x")
        elif self.caseX == IntegrandType.InfiniteInfinite:
            factor = "(Sec[x Pi]^2 Pi / 2) * "
            argX = "Tan[x Pi]"
        # =============== Y ==================
        if self.caseY == IntegrandType.AB:
            factor = factor + "{} * ".format(self.sepY)
            argY = "{} * (y + 1) + {}".format(self.sepY, self.leftY)
        elif self.caseY == IntegrandType.AInfinite:
            factor = factor + self.aInfMapping.FactorString(True, self.leftY, "y")
            argY = self.aInfMapping.ArgString(True, self.leftY, "y")
        elif self.caseY == IntegrandType.InfiniteA:
            factor = factor + self.aInfMapping.FactorString(False, self.rightY, "y")
            argY = self.aInfMapping.ArgString(False, self.rightY, "y")
        elif self.caseY == IntegrandType.InfiniteInfinite:
            factor = factor + "(Sec[y Pi]^2 Pi / 2) * "
            argY = "Tan[y Pi]"
        return "g[x_, y_]:={} f[{}, {}];\n".format(factor, argX, argY)

    def Evaluate(self, x: complex, y: complex) -> complex:
        """
        为了避免大量计算，我们不再使用逼近
        注意和1维不同，2维我们用-1到1的积分
        """
        try:
            newX: complex = x
            newY: complex = y
            factor: complex = 1
            # =============== X ==================
            if self.caseX == IntegrandType.AB:
                factor = factor * self.sepX
                newX = self.sepX * (x + 1) + self.leftX
            elif self.caseX == IntegrandType.AInfinite:
                factor = factor * self.aInfMapping.Factor(True, self.leftX, x)
                newX = self.aInfMapping.Arg(True, self.leftX, x)
            elif self.caseX == IntegrandType.InfiniteA:
                factor = factor * self.aInfMapping.Factor(False, self.rightX, x)
                newX = self.aInfMapping.Arg(False, self.rightX, x)
            elif self.caseX == IntegrandType.InfiniteInfinite:
                b: complex = cmath.pi * x
                factor = factor * 0.5 * cmath.pi / (cmath.cos(b) ** 2)
                newX = cmath.tan(b)
            # =============== Y ==================
            if self.caseY == IntegrandType.AB:
                factor = factor * self.sepY
                newY = self.sepY * (y + 1) + self.leftY
            elif self.caseY == IntegrandType.AInfinite:
                factor = factor * self.aInfMapping.Factor(True, self.leftY, y)
                newY = self.aInfMapping.Arg(True, self.leftY, y)
            elif self.caseY == IntegrandType.InfiniteA:
                factor = factor * self.aInfMapping.Factor(False, self.rightY, y)
                newY = self.aInfMapping.Arg(False, self.rightY, y)
            elif self.caseY == IntegrandType.InfiniteInfinite:
                b: complex = cmath.pi * y
                factor = factor * 0.5 * cmath.pi / (cmath.cos(b) ** 2)
                newY = cmath.tan(b)
            return factor * self.func(newX, newY)
        except (ValueError, ZeroDivisionError):
            return cmath.nan

    def GetDebugInfo(self) -> str:
        head = "Print[\"Original Integral is Integrate[f[x,y], {}x,{},{}{}, {}y,{},{}{}]\"]\n"\
            .format("{", self.leftX, self.rightX, "}", "{", self.leftY, self.rightY, "}")
        if self.hasMathematicaExpression:
            head = head + self.FExpression
        else:
            head = head + self.Dress()
        return head
