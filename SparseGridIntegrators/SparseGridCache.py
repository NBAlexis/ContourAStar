import os

import numpy as np

from SparseGridIntegrators.NestedQuadrature import NestedQuadrature
from SparseGridIntegrators.SparseGridGenerator3D import SparseGrid3D
from SparseGridIntegrators.SparseGridGenerator4D import SparseGrid4D
from SparseGridIntegrators.SparseGridPoints import SparseGridPoints4D, SparseGridPoints3D


def CacheSparseGrid4D(quadrature: NestedQuadrature, folderName: str, maxOrder: int):
    sparseGrid = SparseGrid4D(quadrature)
    [pts, startIndex, endIndex, weights] = sparseGrid.ConstructPointListAndWeightList(maxOrder)
    for i in range(0, maxOrder):
        ptAtOrderI = []
        weightAtOrderI = []
        # ======= get function values ========
        for pointIndex in range(startIndex[i], endIndex[i]):
            ptAtOrderI.append([pts[pointIndex].x, pts[pointIndex].y, pts[pointIndex].z, pts[pointIndex].w])
        for pointIndex in range(0, endIndex[i]):
            weightAtOrderI.append(weights[i][pointIndex])
        np.savetxt(folderName + "/pt" + str(i), ptAtOrderI, delimiter=',')
        np.savetxt(folderName + "/wt" + str(i), weightAtOrderI, delimiter=',')


def CachedPointListAndWeightList(folderName: str, maxOrder: int = 20) -> [list, list, list, list]:
    resPoints = []
    resStartIdx = []
    resEndIdx = []
    resWeights = []
    for i in range(0, maxOrder):
        fileName = folderName + "/pt" + str(i)
        if not os.path.exists(fileName):
            continue
        ptAtOrderI = np.loadtxt(folderName + "/pt" + str(i), delimiter=',').tolist()
        weightAtOrderI = np.loadtxt(folderName + "/wt" + str(i), delimiter=',').tolist()
        if not isinstance(ptAtOrderI[0], list) and 4 == len(ptAtOrderI):
            ptAtOrderI = [ptAtOrderI]
            weightAtOrderI = [weightAtOrderI]
        resStartIdx.append(len(resPoints))
        for pts in ptAtOrderI:
            resPoints.append(SparseGridPoints4D(
                pts[0], pts[1], pts[2], pts[3],
                0, 0, 0, 0, 0, 0, 0, 0))
        resEndIdx.append(len(resPoints))
        resWeights.append(weightAtOrderI)
    return [resPoints, resStartIdx, resEndIdx, resWeights]


def CacheSparseGrid3D(quadrature: NestedQuadrature, folderName: str, maxOrder: int):
    sparseGrid = SparseGrid3D(quadrature)
    [pts, startIndex, endIndex, weights] = sparseGrid.ConstructPointListAndWeightList(maxOrder)
    for i in range(0, maxOrder):
        ptAtOrderI = []
        weightAtOrderI = []
        # ======= get function values ========
        for pointIndex in range(startIndex[i], endIndex[i]):
            ptAtOrderI.append([pts[pointIndex].x, pts[pointIndex].y, pts[pointIndex].z])
        for pointIndex in range(0, endIndex[i]):
            weightAtOrderI.append(weights[i][pointIndex])
        np.savetxt(folderName + "/pt" + str(i), ptAtOrderI, delimiter=',')
        np.savetxt(folderName + "/wt" + str(i), weightAtOrderI, delimiter=',')


def CachedPointListAndWeightList3D(folderName: str, maxOrder: int = 20) -> [list, list, list, list]:
    resPoints = []
    resStartIdx = []
    resEndIdx = []
    resWeights = []
    for i in range(0, maxOrder):
        fileName = folderName + "/pt" + str(i)
        if not os.path.exists(fileName):
            continue
        ptAtOrderI = np.loadtxt(folderName + "/pt" + str(i), delimiter=',').tolist()
        weightAtOrderI = np.loadtxt(folderName + "/wt" + str(i), delimiter=',').tolist()
        if not isinstance(ptAtOrderI[0], list) and 3 == len(ptAtOrderI):
            ptAtOrderI = [ptAtOrderI]
            weightAtOrderI = [weightAtOrderI]
        resStartIdx.append(len(resPoints))
        for pts in ptAtOrderI:
            resPoints.append(SparseGridPoints3D(pts[0], pts[1], pts[2], 0, 0, 0, 0, 0, 0))
        resEndIdx.append(len(resPoints))
        resWeights.append(weightAtOrderI)
    return [resPoints, resStartIdx, resEndIdx, resWeights]
