from UsefulFunctions.CommonFunctions import SaveCSVFile
from UsefulFunctions.NestedQuadrature import NestedQuadrature
from UsefulFunctions.SparseGridGenerator import SparseGrid
from UsefulFunctions.SparseGridGenerator3D import SparseGrid3D


def SaveQuadrature(fileName: str, nestedQuadrature: NestedQuadrature, n: int) -> [list, list]:
    nestedQuadrature.CachePointAndWeights(n)
    [pts3, startIndex3, endIndex3, weights3] = nestedQuadrature.ConstructPointListYLeftMost(n)
    pointAbscissas1D = []
    startEnd1D = []
    for i in range(0, len(pts3)):
        pointAbscissas1D.append([pts3[i].y, pts3[i].yBatch, pts3[i].yIndex])
    for i in range(0, len(startIndex3)):
        startEnd1D.append([startIndex3[i], endIndex3[i]])
    SaveCSVFile("{}D1/Abscissas.csv".format(fileName), pointAbscissas1D, 0, 2)
    SaveCSVFile("{}D1/Index.csv".format(fileName), startEnd1D, 0, 1)
    for i in range(0, len(weights3)):
        weightListSave = [[weights3[i][j]] for j in range(0, len(weights3[i]))]
        SaveCSVFile("{}D1/Weight{}.csv".format(fileName, i), weightListSave, 0, 0)
    sparseGrid2d = SparseGrid(nestedQuadrature)
    [pts2, startIndex2, endIndex2, weights2] = sparseGrid2d.ConstructPointListAndWeightList(n)
    pointAbscissas2D = []
    startEnd2D = []
    for i in range(0, len(pts2)):
        pointAbscissas2D.append([pts2[i].x, pts2[i].y, pts2[i].xBatch, pts2[i].xIndex, pts2[i].yBatch, pts2[i].yIndex])
    for i in range(0, len(startIndex2)):
        startEnd2D.append([startIndex2[i], endIndex2[i]])
    SaveCSVFile("{}D2/Abscissas.csv".format(fileName), pointAbscissas2D, 0, 5)
    SaveCSVFile("{}D2/Index.csv".format(fileName), startEnd2D, 0, 1)
    for i in range(0, len(weights2)):
        weightListSave = [[weights2[i][j]] for j in range(0, len(weights2[i]))]
        SaveCSVFile("{}D2/Weight{}.csv".format(fileName, i), weightListSave, 0, 0)
    sparseGrid = SparseGrid3D(nestedQuadrature)
    [pts, startIndex, endIndex, weights] = sparseGrid.ConstructPointListAndWeightList(n)
    pointAbscissas3D = []
    startEnd3D = []
    for i in range(0, len(pts)):
        pointAbscissas3D.append([pts[i].x, pts[i].y, pts[i].z, pts[i].xBatch, pts[i].xIndex, pts[i].yBatch, pts[i].yIndex, pts[i].zBatch, pts[i].zIndex])
    for i in range(0, len(startIndex)):
        startEnd3D.append([startIndex[i], endIndex[i]])
    SaveCSVFile("{}D3/Abscissas.csv".format(fileName), pointAbscissas3D, 0, 8)
    SaveCSVFile("{}D3/Index.csv".format(fileName), startEnd3D, 0, 1)
    for i in range(0, len(weights)):
        weightListSave = [[weights[i][j]] for j in range(0, len(weights[i]))]
        SaveCSVFile("{}D3/Weight{}.csv".format(fileName, i), weightListSave, 0, 0)

