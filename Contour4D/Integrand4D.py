import cmath

from Contour1D.Integrand import IntegrandType
from Contour2D.Integrand2D import AInfiniteMapping, LogMapping


class Integrand4D:

    def __init__(self, func,
                 leftX: complex, rightX: complex,
                 leftY: complex, rightY: complex,
                 leftZ: complex, rightZ: complex,
                 leftW: complex, rightW: complex,
                 aInfMapping: AInfiniteMapping = LogMapping()):
        self.vpf = None
        self.denominator = None
        self.func = func
        # ================= X ====================
        if cmath.isinf(leftX) and cmath.isinf(rightX):
            self.caseX = IntegrandType.InfiniteInfinite
        elif cmath.isinf(leftX):
            self.caseX = IntegrandType.InfiniteA
        elif cmath.isinf(rightX):
            self.caseX = IntegrandType.AInfinite
        else:
            self.caseX = IntegrandType.AB
        self.leftX: complex = leftX
        self.rightX: complex = rightX
        self.sepX: complex = 0
        if self.caseX == IntegrandType.AB:
            self.sepX = (rightX - leftX) * 0.5
        # ================= Y ====================
        if cmath.isinf(leftY) and cmath.isinf(rightY):
            self.caseY = IntegrandType.InfiniteInfinite
        elif cmath.isinf(leftY):
            self.caseY = IntegrandType.InfiniteA
        elif cmath.isinf(rightY):
            self.caseY = IntegrandType.AInfinite
        else:
            self.caseY = IntegrandType.AB
        self.leftY: complex = leftY
        self.rightY: complex = rightY
        self.sepY: complex = 0
        if self.caseY == IntegrandType.AB:
            self.sepY = (rightY - leftY) * 0.5
        # ================= Z ====================
        if cmath.isinf(leftZ) and cmath.isinf(rightZ):
            self.caseZ = IntegrandType.InfiniteInfinite
        elif cmath.isinf(leftZ):
            self.caseZ = IntegrandType.InfiniteA
        elif cmath.isinf(rightZ):
            self.caseZ = IntegrandType.AInfinite
        else:
            self.caseZ = IntegrandType.AB
        self.leftZ: complex = leftZ
        self.rightZ: complex = rightZ
        self.sepZ: complex = 0
        if self.caseZ == IntegrandType.AB:
            self.sepZ = (rightZ - leftZ) * 0.5
        # ================= W ====================
        if cmath.isinf(leftW) and cmath.isinf(rightW):
            self.caseW = IntegrandType.InfiniteInfinite
        elif cmath.isinf(leftW):
            self.caseW = IntegrandType.InfiniteA
        elif cmath.isinf(rightW):
            self.caseW = IntegrandType.AInfinite
        else:
            self.caseW = IntegrandType.AB
        self.leftW: complex = leftW
        self.rightW: complex = rightW
        self.sepW: complex = 0
        if self.caseW == IntegrandType.AB:
            self.sepW = (rightW - leftW) * 0.5
        # ================ Others ==================
        self.aInfMapping = aInfMapping
        self.hasMathematicaExpression = False
        self.FExpression = ""

    def MakeSureZeroOne(self):
        self.caseX = IntegrandType.ZeroOne
        self.caseY = IntegrandType.ZeroOne
        self.caseZ = IntegrandType.ZeroOne
        self.caseW = IntegrandType.ZeroOne

    def Evaluate(self, x: complex, y: complex, z: complex, w: complex) -> complex:
        """
        为了避免大量计算，我们不再使用逼近
        注意和1维不同，2维我们用-1到1的积分
        """
        try:
            newX: complex = x
            newY: complex = y
            newZ: complex = z
            newW: complex = w
            factor: complex = 1
            # =============== X ==================
            if self.caseX == IntegrandType.AB or self.caseX == IntegrandType.ZeroOne:
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
            if self.caseY == IntegrandType.AB or self.caseY == IntegrandType.ZeroOne:
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
            # =============== Z ==================
            if self.caseZ == IntegrandType.AB or self.caseZ == IntegrandType.ZeroOne:
                factor = factor * self.sepZ
                newZ = self.sepZ * (z + 1) + self.leftZ
            elif self.caseZ == IntegrandType.AInfinite:
                factor = factor * self.aInfMapping.Factor(True, self.leftZ, z)
                newZ = self.aInfMapping.Arg(True, self.leftZ, z)
            elif self.caseZ == IntegrandType.InfiniteA:
                factor = factor * self.aInfMapping.Factor(False, self.rightZ, z)
                newZ = self.aInfMapping.Arg(False, self.rightZ, z)
            elif self.caseZ == IntegrandType.InfiniteInfinite:
                b: complex = cmath.pi * z
                factor = factor * 0.5 * cmath.pi / (cmath.cos(b) ** 2)
                newZ = cmath.tan(b)
            # =============== W ==================
            if self.caseW == IntegrandType.AB or self.caseW == IntegrandType.ZeroOne:
                factor = factor * self.sepW
                newW = self.sepW * (w + 1) + self.leftW
            elif self.caseW == IntegrandType.AInfinite:
                factor = factor * self.aInfMapping.Factor(True, self.leftW, w)
                newW = self.aInfMapping.Arg(True, self.leftW, w)
            elif self.caseW == IntegrandType.InfiniteA:
                factor = factor * self.aInfMapping.Factor(False, self.rightW, w)
                newW = self.aInfMapping.Arg(False, self.rightW, w)
            elif self.caseW == IntegrandType.InfiniteInfinite:
                b: complex = cmath.pi * w
                factor = factor * 0.5 * cmath.pi / (cmath.cos(b) ** 2)
                newW = cmath.tan(b)
            return factor * self.func(newX, newY, newZ, newW)
        except (ValueError, ZeroDivisionError):
            return cmath.nan

    def Dress(self) -> str:
        argX = ""
        argY = ""
        argZ = ""
        argW = ""
        factor = ""
        # =============== X ==================
        if self.caseX == IntegrandType.ZeroOne:
            factor = "(1/2) * "
            argX = "(x + 1)/2"
        elif self.caseX == IntegrandType.AB:
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
        if self.caseY == IntegrandType.ZeroOne:
            factor = factor + "(1/2) * "
            argY = "(y + 1)/2"
        elif self.caseY == IntegrandType.AB:
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
        # =============== Z ==================
        if self.caseZ == IntegrandType.ZeroOne:
            factor = factor + "(1/2) * "
            argZ = "(z + 1)/2"
        elif self.caseZ == IntegrandType.AB:
            factor = factor + "{} * ".format(self.sepZ)
            argZ = "{} * (z + 1) + {}".format(self.sepZ, self.leftZ)
        elif self.caseZ == IntegrandType.AInfinite:
            factor = factor + self.aInfMapping.FactorString(True, self.leftZ, "z")
            argZ = self.aInfMapping.ArgString(True, self.leftZ, "z")
        elif self.caseZ == IntegrandType.InfiniteA:
            factor = factor + self.aInfMapping.FactorString(False, self.rightZ, "z")
            argZ = self.aInfMapping.ArgString(False, self.rightZ, "z")
        elif self.caseZ == IntegrandType.InfiniteInfinite:
            factor = factor + "(Sec[z Pi]^2 Pi / 2) * "
            argZ = "Tan[z Pi]"
        # =============== W ==================
        if self.caseW == IntegrandType.ZeroOne:
            factor = factor + "(1/2) * "
            argW = "(w + 1)/2"
        elif self.caseW == IntegrandType.AB:
            factor = factor + "{} * ".format(self.sepW)
            argW = "{} * (w + 1) + {}".format(self.sepW, self.leftW)
        elif self.caseW == IntegrandType.AInfinite:
            factor = factor + self.aInfMapping.FactorString(True, self.leftW, "w")
            argW = self.aInfMapping.ArgString(True, self.leftW, "w")
        elif self.caseW == IntegrandType.InfiniteA:
            factor = factor + self.aInfMapping.FactorString(False, self.rightW, "w")
            argW = self.aInfMapping.ArgString(False, self.rightW, "w")
        elif self.caseW == IntegrandType.InfiniteInfinite:
            factor = factor + "(Sec[w Pi]^2 Pi / 2) * "
            argZ = "Tan[w Pi]"
        return "g[x_, y_, z_, w_]:={} f[{}, {}, {}, {}];".format(factor, argX, argY, argZ, argW)

    def GetDebugInfo(self) -> str:
        head = "Print[\"Original Integral is Integrate[f[x,y,z,w], {}x,{},{}{}, {}y,{},{}{}, {}z,{},{}{}, {}w,{},{}{}]\"]\n" \
            .format("{", self.leftX, self.rightX, "}",
                    "{", self.leftY, self.rightY, "}",
                    "{", self.leftZ, self.rightZ, "}",
                    "{", self.leftW, self.rightW, "}")
        if self.hasMathematicaExpression:
            return head + self.FExpression
        return head + self.Dress()

    def SetMathematicaExpress(self, expression: str) -> str:
        self.hasMathematicaExpression = True
        self.FExpression = "f[x_, y_, z_, w_]:=" + expression + ";\n" + self.Dress() + "\n"
        return self.FExpression

    def GetMathematicaExpress(self) -> [bool, str]:
        if self.hasMathematicaExpression:
            return [True, self.FExpression]
        return [False, ""]

    def SetVPF(self, vpf):
        self.vpf = vpf

    def GetVPF(self):
        return self.vpf

    def SetDenominator(self, denominator):
        self.denominator = denominator

    def HasDenominator(self) -> bool:
        return self.denominator is not None

    def EvaluateDenominator(self, x: complex, y: complex, z: complex, w: complex) -> complex:
        """
        it should be a polynomial, so we do not use "try"
        """
        newX: complex = x
        newY: complex = y
        newZ: complex = z
        newW: complex = w
        # =============== X ==================
        if self.caseX == IntegrandType.AB or self.caseX == IntegrandType.ZeroOne:
            newX = self.sepX * (x + 1) + self.leftX
        elif self.caseX == IntegrandType.AInfinite:
            newX = self.aInfMapping.Arg(True, self.leftX, x)
        elif self.caseX == IntegrandType.InfiniteA:
            newX = self.aInfMapping.Arg(False, self.rightX, x)
        elif self.caseX == IntegrandType.InfiniteInfinite:
            b: complex = cmath.pi * x
            newX = cmath.tan(b)
        # =============== Y ==================
        if self.caseY == IntegrandType.AB or self.caseY == IntegrandType.ZeroOne:
            newY = self.sepY * (y + 1) + self.leftY
        elif self.caseY == IntegrandType.AInfinite:
            newY = self.aInfMapping.Arg(True, self.leftY, y)
        elif self.caseY == IntegrandType.InfiniteA:
            newY = self.aInfMapping.Arg(False, self.rightY, y)
        elif self.caseY == IntegrandType.InfiniteInfinite:
            b: complex = cmath.pi * y
            newY = cmath.tan(b)
        # =============== Z ==================
        if self.caseZ == IntegrandType.AB or self.caseZ == IntegrandType.ZeroOne:
            newZ = self.sepZ * (z + 1) + self.leftZ
        elif self.caseZ == IntegrandType.AInfinite:
            newZ = self.aInfMapping.Arg(True, self.leftZ, z)
        elif self.caseZ == IntegrandType.InfiniteA:
            newZ = self.aInfMapping.Arg(False, self.rightZ, z)
        elif self.caseZ == IntegrandType.InfiniteInfinite:
            b: complex = cmath.pi * z
            newZ = cmath.tan(b)
        # =============== W ==================
        if self.caseW == IntegrandType.AB or self.caseW == IntegrandType.ZeroOne:
            newW = self.sepW * (w + 1) + self.leftW
        elif self.caseW == IntegrandType.AInfinite:
            newW = self.aInfMapping.Arg(True, self.leftW, w)
        elif self.caseW == IntegrandType.InfiniteA:
            newW = self.aInfMapping.Arg(False, self.rightW, w)
        elif self.caseW == IntegrandType.InfiniteInfinite:
            b: complex = cmath.pi * w
            newW = cmath.tan(b)
        return self.denominator(newX, newY, newZ, newW)