"""
    Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
    Last update: 02/12/2014
""" 
from _hamiltonian import HamiltonianParams
from _hamiltonian import Hamiltonian
from _hamiltonian import CxHamiltonian

from _hamiltonian import GrapheneTbParams
from _hamiltonian import GrapheneTbHam

from _hamiltonian import GrapheneKpParams
from _hamiltonian import GrapheneKpHam

from _hamiltonian import TISurfKpParams
from _hamiltonian import TISurfKpHam

from _hamiltonian import TISurfKpParams4
from _hamiltonian import TISurfKpHam4

from qmicad.atoms import PeriodicTable

# TI surface k.p default parameters
_TISurfKpParamsOrgInit =  TISurfKpParams.__init__
def _TISurfKpParamsInit(self):
    _TISurfKpParamsOrgInit(self)
    self.dtol   = 1E-3              # Tolerance when considering neighbors
    self.a      = 2.0               # Lattice spacing
    self.K      = 1.165             # Coefficient ot solve Fermion doubling
    self.A2     = 3.33              # A2 paramter
    self.C      = 0.0               # C parameter
    self.ptable = PeriodicTable()   # Periodic table for k.p for TI
    self.ptable.add(0, "D", 2, 2)   # Fake atom for TI k.p
    
TISurfKpParams.__init__ = _TISurfKpParamsInit

# TI surface k.p pickle support
def _TISurfKpParamsSetState(self, dct):
    self.__dict__.update(dct)
#    self.update()
    
TISurfKpParams.__setstate__ = _TISurfKpParamsSetState

def _TISurfKpParamsGetState(self):
    dct = dict(self.__dict__)
    return dct

TISurfKpParams.__getstate__ = _TISurfKpParamsGetState
TISurfKpParams.__getstate_manages_dict__ = True

# TI surface k.p default parameters in four spins basis set.
_TISurfKpParams4OrgInit =  TISurfKpParams4.__init__
def _TISurfKpParams4Init(self):
    _TISurfKpParams4OrgInit(self)
    self.dtol   = 1E-3              # Tolerance when considering neighbors
    self.a      = 2.0               # Lattice spacing
    self.K      = 1.165             # Coefficient ot solve Fermion doubling
    self.A2     = 3.33              # A2 paramter
    self.C      = 0.0               # C parameter
    self.ptable = PeriodicTable()   # Periodic table for k.p for TI
    self.ptable.add(0, "D", 4, 4)   # Fake atom for TI k.p
    
TISurfKpParams4.__init__ = _TISurfKpParams4Init

# TI surface k.p pickle support
def _TISurfKpParams4SetState(self, dct):
    self.__dict__.update(dct)
#    self.update()
    
TISurfKpParams4.__setstate__ = _TISurfKpParams4SetState

def _TISurfKpParams4GetState(self):
    dct = dict(self.__dict__)
    return dct

TISurfKpParams4.__getstate__ = _TISurfKpParams4GetState
TISurfKpParams4.__getstate_manages_dict__ = True

# Graphene k.p default parameters
_GrapheneKpParamsOrgInit =  GrapheneKpParams.__init__
def _GrapheneKpParamsInit(self):
    _GrapheneKpParamsOrgInit(self)
    self.dtol   = 1E-3              # Tolerance when considering neighbors
    self.a      = 4.0               # Lattice spacing
    self.K      = 1.165             # Coefficient ot solve Fermion doubling
    self.gamma  = 3.16*1.42*3/2     # gamma = hbar * v_F
    self.ptable = PeriodicTable()   # Periodic table for graphene k.p
    self.ptable.add(0, "D", 2, 2)   # Fake atom for k.p
    
GrapheneKpParams.__init__ = _GrapheneKpParamsInit

# Graphene k.p pickle support
def _GrapheneKpParamsSetState(self, dct):
    self.__dict__.update(dct)
#    self.update()
    
GrapheneKpParams.__setstate__ = _GrapheneKpParamsSetState

def _GrapheneKpParamsGetState(self):
    dct = dict(self.__dict__)
    return dct

GrapheneKpParams.__getstate__ = _GrapheneKpParamsGetState
GrapheneKpParams.__getstate_manages_dict__ = True


# Graphene tight binding default parameters
_GrapheneTbParamsOrgInit =  GrapheneTbParams.__init__
def _GrapheneTbParamsInit(self):
    _GrapheneTbParamsOrgInit(self)
    self.dtol   = 1E-3              # Tolerance when considering neighbors
    self.ec         = 0             # Site energy of C
    self.di0        = 1.42          # In-plane C-C bond length 
    self.ti0        = 3.16          # In plane C-C hopping parameter
    self.do0        = 3.35          # Out-of-plane C-C distance
    self.to0        = 0.39          # Out-of-pane C-C hopping parameter
    self.lmdz       = 0.6           # For interlayer all nearest neighbors
    self.lmdxy      = 1.7           # See PRL 109, 236604 (2012)
    self.alpha      = 1.65
    self.doX        = 6.0           # Out-of-plane neighbor cut-off distance
    self.ptable = PeriodicTable()   # Periodic table for graphene
    
GrapheneTbParams.__init__ = _GrapheneTbParamsInit

# Graphene tight binding pickle support
def _GrapheneTbParamsSetState(self, dct):
    self.__dict__.update(dct)
#    self.update()
    
GrapheneTbParams.__setstate__ = _GrapheneTbParamsSetState

def _GrapheneTbParamsGetState(self):
    dct = dict(self.__dict__)
    return dct

GrapheneTbParams.__getstate__ = _GrapheneTbParamsGetState
GrapheneTbParams.__getstate_manages_dict__ = True


