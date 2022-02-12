from SparseGridIntegrators.ClenshawCurtis import ClenshawCurtis
from SparseGridIntegrators.ClenshawCurtisExp import ClenshawCurtisExp
from SparseGridIntegrators.GaussianPatterson import GaussPatterson
from SparseGridIntegrators.NestedQuadrature import Trapezoidal
from SparseGridIntegrators.SparseGridCache import CacheSparseGrid4D, CacheSparseGrid3D

CacheSparseGrid4D(GaussPatterson(), "../_Data/SparseGrid4D/GaussPatterson", 9)
CacheSparseGrid4D(ClenshawCurtis(), "../_Data/SparseGrid4D/ClenshawCurtis", 18)
CacheSparseGrid4D(ClenshawCurtisExp(), "../_Data/SparseGrid4D/ClenshawCurtisExp", 12)
CacheSparseGrid4D(Trapezoidal(), "../_Data/SparseGrid4D/Trapezoidal", 7)

CacheSparseGrid3D(GaussPatterson(), "../_Data/SparseGrid3D/GaussPatterson", 9)
CacheSparseGrid3D(ClenshawCurtis(), "../_Data/SparseGrid3D/ClenshawCurtis", 18)
CacheSparseGrid3D(ClenshawCurtisExp(), "../_Data/SparseGrid3D/ClenshawCurtisExp", 12)
CacheSparseGrid3D(Trapezoidal(), "../_Data/SparseGrid3D/Trapezoidal", 7)
