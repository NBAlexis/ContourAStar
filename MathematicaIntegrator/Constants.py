import math

from wolframclient.language import wl, wlexpr
from wolframclient.evaluation import WolframLanguageSession
from wolframclient.deserializers import WXFConsumer

_MATHEMATICA_PATH_ = 'F:\\Mathematica\\12.0\\MathKernel.exe'

Complex = wl.Complex


class MathConsumer(WXFConsumer):
    """Implement a consumer with basic arithmetic operation."""

    # Specific convertion for Pi, other symbols use the default method.
    def consume_symbol(self, current_token, tokens, **kwargs):
        # Convert symbol Pi to its numeric value as defined in Python
        if current_token.data == 'Pi':
            return math.pi
        else:
            return super().consume_symbol(current_token, tokens, **kwargs)

    # Associate heads with the method to convert them to Python types.
    DISPATCH = {
        Complex: 'build_complex',
        wl.Rational: 'build_rational',
        wl.Plus: 'build_plus',
        wl.Times: 'build_times'
    }

    # Overload the method that builds functions.
    def build_function(self, head, args, **kwargs):
        # check if there is a specific function associated to the function head
        builder_func = self.DISPATCH.get(head, None)
        if builder_func is not None:
            try:
                # get the class method and apply it to the arguments.
                return getattr(self, builder_func)(*args)
            except Exception:
                # instead of failing, fallback to default case.
                return super().build_function(head, args, **kwargs)
        # heads not listed in DISPATCH are delegated to parent's method
        else:
            return super().build_function(head, args, **kwargs)

    def build_plus(self, *args):
        total = 0
        for arg in args:
            total = total + arg
        return total

    def build_times(self, *args):
        total = 1
        for arg in args:
            total = total * arg
        return total

    def build_rational(self, *args):
        if len(args) != 2:
            raise ValueError('Rational format not supported.')
        return args[0] / args[1]

    def build_complex(self, *args):
        if len(args) != 2:
            raise ValueError('Complex format not supported.')
        return complex(args[0], args[1])


class MathLink:

    def __init__(self):
        self.session = WolframLanguageSession(_MATHEMATICA_PATH_)

    def Call(self, cmd: str):
        return self.session.evaluate(wlexpr(cmd))

    def Quit(self):
        self.session.stop()


def CallMath():
    session = WolframLanguageSession(_MATHEMATICA_PATH_)
    [_, resr, resi] = session.evaluate(wlexpr('Quiet[If[res = Check[NIntegrate[1/(x + 0.5 I), {x, 0, 1}], False, {NIntegrate::slwcon, NIntegrate::ncvb}]; BooleanQ[res], {False, 0, 0}, {True, Re[res], Im[res]}], {NIntegrate::slwcon, NIntegrate::ncvb}]'))
    # complex_result = binary_deserialize(res, consumer=MathConsumer())
    print(resr)
    print(resi)
    session.stop()
