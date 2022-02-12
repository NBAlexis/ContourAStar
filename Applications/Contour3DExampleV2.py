import numpy as np

from Contour1D.CommonDefinitions import LogLevel
from Contour1D.IntegrandV2 import IntegrandV2
from Contour3D.ExtendPath3DV2 import ExtendPath3DV2
from Contour3D.Integrator3DV2 import SparseGridIntegrator3DV2


def intfunc1(x: complex, y: complex, z: complex) -> complex:
    return 1 / (np.log((x + y + z)) + 0.5)


def intfunc1d(x: complex, y: complex, z: complex) -> complex:
    return np.log((x + y + z)) + 0.5


intgrand1 = IntegrandV2(intfunc1, intfunc1d, "1 / (Log[x + y + z] + 0.5)", "f[x,y,z]")


def intfunc2(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 / ((x + y + z - 0.25) ** 2)


def intfunc2d(x: complex, y: complex, z: complex) -> complex:
    return 2.5 + x + y + z


intgrand2 = IntegrandV2(intfunc2, intfunc2d, "0.5/(2.5 + x + y + z)^2", "f[x,y,z]")


def intfunc3(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 / ((x + y + z + x * x - 0.25) ** 2)


def intfunc3d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return x + y + z + x * x - 0.25


intgrand3 = IntegrandV2(intfunc3, intfunc3d, "1 / ((x + y + z + x * x - 0.25) ^ 2)",
                        "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc4(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 / ((0.1 * x + 0.2 * y + 0.3 * z
                     + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                     + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                     + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ** 2)


def intfunc4d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return (0.1 * x + 0.2 * y + 0.3 * z
            + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
            + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
            + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5)


intgrand4 = IntegrandV2(intfunc4, intfunc4d, """1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ^ 2)""",
                        "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc5(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 / ((0.1 * x + 0.2 * y + 0.3 * z
                     + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                     + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                     + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                     - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                     + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                     + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                     + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                     + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2)


def intfunc5d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return (0.1 * x + 0.2 * y + 0.3 * z
            + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
            + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
            + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
            - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
            + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
            + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
            + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
            + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x)


intgrand5 = IntegrandV2(intfunc5, intfunc5d, """1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2)""",
                        "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc6(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 / ((0.1 * x + 0.2 * y + 0.3 * z
                     + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                     + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                     + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                     - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                     + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                     + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                     + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                     - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2)


def intfunc6d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return (0.1 * x + 0.2 * y + 0.3 * z
            + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
            + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
            + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
            - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
            + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
            + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
            + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
            - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x)


intgrand6 = IntegrandV2(intfunc6, intfunc6d, """1 / ((0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2)""",
                        "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc7(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 * (1 - x) * (1 - x) * (1 - y) / ((-0.02 - (1 - x) * (1 - y) * x * z
                                                   + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)


def intfunc7d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return -0.02 - (1 - x) * (1 - y) * x * z + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)


intgrand7 = IntegrandV2(intfunc7, intfunc7d, """(1 - x) * (1 - x) * (1 - y) / ((-0.02 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ^ 2)""",
                        "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc8(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 * (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                                   + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)


def intfunc8d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.2 - (1 - x) * (1 - y) * x * z + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)


intgrand8 = IntegrandV2(intfunc8, intfunc8d, """(1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ^ 2)""",
                        "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc9(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 * (1 - x) * (1 - x) * (1 - y) / ((0.2 - 0.142857 * x * (1 - x) * y
                                                   - x * (1 - x) * (1 - y) * z
                                                   - 0.142857 * (1 - x) * (1 - x) * y * (1 - y) * z
                                                   + 0.4 * (1 - x) * (1 - x) * y * (1 - y) * (1 - z)) ** 2)


def intfunc9d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return (0.2 - 0.142857 * x * (1 - x) * y
            - x * (1 - x) * (1 - y) * z
            - 0.142857 * (1 - x) * (1 - x) * y * (1 - y) * z
            + 0.4 * (1 - x) * (1 - x) * y * (1 - y) * (1 - z))


intgrand9 = IntegrandV2(intfunc9, intfunc9d, """(1-x)*(1-x)*(1-y)/((0.2 - 0.142857*x*(1-x)*y
                               - x*(1-x)*(1-y)*z
                               - 0.142857*(1-x)*(1-x)*y*(1-y)*z 
                               + 0.4*(1-x)*(1-x)*y*(1-y)*(1-z))^2)""",
                        "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc10(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 * (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                                   - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ** 2)


def intfunc10d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.2 - (1 - x) * (1 - y) * x * z - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y


intgrand10 = IntegrandV2(intfunc10, intfunc10d, """(1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ^ 2)""",
                         "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc11(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 * (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                                   + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ** 2)


def intfunc11d(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.2 - (1 - x) * (1 - y) * x * z + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y


intgrand11 = IntegrandV2(intfunc11, intfunc11d, """(1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                           + 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z) * y) ^ 2)""",
                         "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")


def intfunc12(x: complex, y: complex, z: complex) -> complex:
    x = 0.5 * (x + 1)
    y = 0.5 * (y + 1)
    z = 0.5 * (z + 1)
    return 0.125 * (1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                                   - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ** 2)


# def intfunc12d(x: complex, y: complex, z: complex) -> complex:
#     x = 0.5 * (x + 1)
#     y = 0.5 * (y + 1)
#     z = 0.5 * (z + 1)
#     return 0.2 - (1 - x) * (1 - y) * x * z - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)


intgrand12 = IntegrandV2(intfunc12, None, """(1 - x) * (1 - x) * (1 - y) / ((0.2 - (1 - x) * (1 - y) * x * z
                                                   - 0.4 * (1 - x) * (1 - x) * (1 - y) * (1 - z)) ^ 2)""",
                         "f[(x+1)/2,(y+1)/2,(z+1)/2]/8")

intgrator = SparseGridIntegrator3DV2(logLevel=LogLevel.Warning)
grid3d = ExtendPath3DV2(5, 5, 0, intgrator, logLevel=LogLevel.Warning)
grid3dmid = ExtendPath3DV2(7, 7, 0, intgrator, logLevel=LogLevel.Warning)
grid3d9 = ExtendPath3DV2(9, 9, 0, intgrator, logLevel=LogLevel.Warning)

res = ""

print(grid3d.Integrate(intgrand1))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand2))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand3))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand4))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand5))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand6))
res = res + grid3d.GatherInfo()

print(grid3dmid.Integrate(intgrand7))
res = res + grid3dmid.GatherInfo()

print(grid3d.Integrate(intgrand8))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand9))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand10))
res = res + grid3d.GatherInfo()

print(grid3d.Integrate(intgrand11))
res = res + grid3d.GatherInfo()

print(grid3d9.Integrate(intgrand12))
res = res + grid3d9.GatherInfo()

print(res)

"""
[<AStarResult.Finished: 1>, (1.7395625773125973+6.186906021257622j)]
[<AStarResult.Finished: 1>, (1.1499386328424233-0.7853981633972124j)]
[<AStarResult.Finished: 1>, (0.8447069460960448-0.6506451422834971j)]
[<AStarResult.Finished: 1>, (-2.849446283400229+0.36359972372883664j)]
[<AStarResult.Finished: 1>, (-3.4488471218798096-1.6557182004182118j)]
[<AStarResult.Finished: 1>, (-3.6329746931687907+1.7539244081342726j)]
[<AStarResult.Finished: 1>, (-19.495544330786263-17.829300247607733j)]
[<AStarResult.Finished: 1>, (4.64119089459386-4.5105233894253995j)]
[<AStarResult.Finished: 1>, (7.120704446607922-8.8990564086542j)]
[<AStarResult.Finished: 1>, (10.05697561629099-7.534126780833717j)]
[<AStarResult.Finished: 1>, (6.5237277112999426-6.591009777256285j)]
[<AStarResult.Finished: 1>, (16.74584627545428+3.5342681046993567j)]
"""
