"""  
 Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 Last update: 02/12/2014
"""

from _utils import Printable
from _utils import Option
from _utils import Timer
from _utils import CxMatp
from _utils import Matp
from _utils import Vecp
from _utils import VecGrid
from _utils import Workers
from _utils import Point
from _utils import Quadrilateral

##
# Enum pickler.
#
def isEnumType(o):
    return isinstance(o, type) and issubclass(o,int) and not (o is int)

def _tuple2enum(enum, value):
    enum = getattr(_utils, enum)
    e = enum.values.get(value,None)
    if e is None:
        e = enum(value)
    return e

def _registerEnumPicklers(): 
    from copy_reg import constructor, pickle
    def reduce_enum(e):
        enum = type(e).__name__.split('.')[-1]
        return ( _tuple2enum, ( enum, int(e) ) )
    constructor( _tuple2enum)
    for e in [ e for e in vars(_utils).itervalues() if isEnumType(e) ]:
        pickle(e, reduce_enum)

_registerEnumPicklers()
