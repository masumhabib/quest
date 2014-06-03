""" 
    Non-equilibrium Green Function Formalism simulation package. It contains the 
    following:
        * Coherent Recursive Green Function Algorithm (RGFA)
        * Parallel Energy Loop
    
    Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
    Last update: 02/12/2014
"""
 
from _negf import CohRgfaParams
from _negf import NegfEloop

# Coherent RGFA default parameters
_CohRgfaParamsOrgInit =  CohRgfaParams.__init__
def _CohRgfaParamsInit(self, nb):
    _CohRgfaParamsOrgInit(self, nb)
    self.kT             = 0.0259        # Temperature in eV (300 K)
    self.ieta           = 1E-3j         # Contact imaginary potential
    self.mu             = 0.0           # Device Fermi level
    self.muS            = 0.0           # Source Fermi level
    self.muD            = 0.0           # Source Fermi level
    self.isOrthogonal   = True          # Orthogonal basis?
    self.Emin           =-1.0           # Minimum energy 
    self.Emax           = 1.0           # Maximum energy
    self.dE             = 0.005         # Energy step
    self.AutoGenE       = False         # Generate grid automatically?
    self.AdaptiveGrid   = False         # Adaptive E grid?
    
CohRgfaParams.__init__ = _CohRgfaParamsInit

_CohRgfaParamsOrgStr =  CohRgfaParams.__str__
def _CohRgfaParamsStr(self):
    msg = _CohRgfaParamsOrgStr(self)
    msg += "\n" + self.Prefix + " Emin         = " + str(self.Emin)
    msg += "\n" + self.Prefix + " Emax         = " + str(self.Emax)
    msg += "\n" + self.Prefix + " dE           = " + str(self.dE)
    msg += "\n" + self.Prefix + " AutoGenE     = " + str(self.AutoGenE)
    msg += "\n" + self.Prefix + " AdaptiveGrid = " + str(self.AdaptiveGrid)
    
    return msg
CohRgfaParams.__str__ = _CohRgfaParamsStr

# Pickle support
def _CohRgfaParamsSetState(self, dct):
    self.__dict__.update(dct)
    
CohRgfaParams.__setstate__ = _CohRgfaParamsSetState

def _CohRgfaParamsGetState(self):
    dct = dict(self.__dict__)
    return dct

CohRgfaParams.__getstate__ = _CohRgfaParamsGetState
CohRgfaParams.__getstate_manages_dict__ = True


def _CohRgfaParamsGetInitArgs(self):
    return (self.nb,)

CohRgfaParams.__getinitargs__ = _CohRgfaParamsGetInitArgs
