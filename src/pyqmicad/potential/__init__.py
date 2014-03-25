##  
#  Author: K M Masum Habib <masum.habib@virginia.edu>
#  Last update: 02/12/2014
 
from _potential import Potential
from _potential import LinearPot

# Electrostatic potential
_LinearPotOrgInit =  LinearPot.__init__
def _LinearPotInit(self, atoms = None):
    _LinearPotOrgInit(self, atoms)
    self.rVS        =-0.5          # Fraction of voltage to be applied to drain
    self.rVD        = 0.5          # Fraction of voltage to be applied to drain
    self.rVG        = []           # Gate voltage ratios
    self.gql        = []           # Gate quadrilaterals just for pickle
    self.lql        = []           # Linear region quadrilaterals just for pickle

LinearPot.__init__ = _LinearPotInit

# Pickle support
def _LinearPotSetState(self, dct):
    self.__dict__.update(dct)
    
LinearPot.__setstate__ = _LinearPotSetState

def _LinearPotGetState(self):
    dct = dict(self.__dict__)
    return dct

LinearPot.__getstate__ = _LinearPotGetState
LinearPot.__getstate_manages_dict__ = True

