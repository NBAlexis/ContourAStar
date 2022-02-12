import numpy as np


class IntegrandV2:

    def __init__(self, func, denom, sFunc, sIntegrand):
        self.vFunc = np.vectorize(func)
        self.vDenom = None if denom is None else np.vectorize(denom)
        self.sFunc = sFunc
        self.sIntegrand = sIntegrand
