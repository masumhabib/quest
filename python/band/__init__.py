"""
    Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
    Last update: 02/12/2014
"""

from _band import BandStructParams
from _band import BandStruct

# Coherent RGFA default parameters
_BandStructParamsOrgInit =  BandStructParams.__init__
def _BandStructParamsInit(self, nn):
    _BandStructParamsOrgInit(self, nn)
    self.isOrthogonal   = True          # Orthogonal basis?
    
BandStructParams.__init__ = _BandStructParamsInit