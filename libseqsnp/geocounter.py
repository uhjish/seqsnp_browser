import sys
sys.path.append('/usr/local/lib')

from pygeoseqlib.query import GeoseqQuery
from pygeoseqlib.count import PyGeocount


class GeoTileCounter:
    def __init__(self,libfile):
        self.libfile = libfile
        self.program = "count"
        self.gq = GeoseqQuery(self.program,"",0,0,self.libfile)
    def getcount(self,seq):
        ts = len(seq)
        self.gq.tile_size=ts
        ret=self.gq.execute_single_query(seq)
        fwd = ret.total_forward
        rev = ret.total_reverse
        return (fwd+rev)

