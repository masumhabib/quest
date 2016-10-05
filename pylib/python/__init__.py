"""
 QUantum mechanics Enabled Simulation Toolset (QUEST) library is a 
 collection of C++ classes and functions for simulation and design of 
 nano-scaled devices using quantum mechanical tools such as Non-equilibrium 
 Green's Function (NEGF) formalism. This library provides a common framework 
 for NEGF so that it can be used with any empirical tight binding, k.p model, 
 extended Huckel method and density functional theory codes. This library also
 includes a generic empirical tight binding model and a k.p model which can be 
 extended for any material with known parameters.
  
 Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 Last update: 02/12/2014

"""

import DLFCN as _dl
import sys as _sys
if _sys.platform == "linux" or _sys.platform == "linux2":
    _sys.setdlopenflags(_dl.RTLD_NOW | _dl.RTLD_GLOBAL)

from quest import * 
import tmfsc
import linspace 
import vprint
import simulators

# TI surface k.p default parameters
_TISurfKpParamsOrgInit =  hamiltonian.TISurfKpParams.__init__
def _TISurfKpParamsInit(self):
    _TISurfKpParamsOrgInit(self)
#    self.dtol   = 1E-3              # Tolerance when considering neighbors
#    self.a      = 2.0               # Lattice spacing
#    self.K      = 1.0               # Coefficient ot solve Fermion doubling
#    self.A2     = 3.33              # A2 paramter
#    self.C      = 0.0               # C parameter
#    self.ptable = atoms.PeriodicTable()   # Periodic table for k.p for TI
#    self.ptable.add(0, "D", 2, 2)   # Fake atom for TI k.p
    
hamiltonian.TISurfKpParams.__init__ = _TISurfKpParamsInit


# TI surface k.p pickle support
def _TISurfKpParamsSetState(self, dct):
    self.__dict__.update(dct)
    
hamiltonian.TISurfKpParams.__setstate__ = _TISurfKpParamsSetState

def _TISurfKpParamsGetState(self):
    dct = dict(self.__dict__)
    return dct

hamiltonian.TISurfKpParams.__getstate__ = _TISurfKpParamsGetState
hamiltonian.TISurfKpParams.__getstate_manages_dict__ = True

# TI surface k.p default parameters in four spins basis set.
_TISurfKpParams4OrgInit =  hamiltonian.TISurfKpParams4.__init__
def _TISurfKpParams4Init(self):
    _TISurfKpParams4OrgInit(self)
#    self.dtol   = 1E-3              # Tolerance when considering neighbors
#    self.a      = 2.0               # Lattice spacing
#    self.K      = 1.165             # Coefficient ot solve Fermion doubling
#    self.A2     = 3.33              # A2 paramter
#    self.C      = 0.0               # C parameter
#    self.ptable = atoms.PeriodicTable()   # Periodic table for k.p for TI
#    self.ptable.add(0, "D", 4, 4)   # Fake atom for TI k.p
    
hamiltonian.TISurfKpParams4.__init__ = _TISurfKpParams4Init

# TI surface k.p pickle support
def _TISurfKpParams4SetState(self, dct):
    self.__dict__.update(dct)
    
hamiltonian.TISurfKpParams4.__setstate__ = _TISurfKpParams4SetState

def _TISurfKpParams4GetState(self):
    dct = dict(self.__dict__)
    return dct

hamiltonian.TISurfKpParams4.__getstate__ = _TISurfKpParams4GetState
hamiltonian.TISurfKpParams4.__getstate_manages_dict__ = True

# Graphene k.p default parameters
_GrapheneKpParamsOrgInit =  hamiltonian.GrapheneKpParams.__init__
def _GrapheneKpParamsInit(self):
    _GrapheneKpParamsOrgInit(self)
#    self.dtol   = 1E-3              # Tolerance when considering neighbors
#    self.a      = 4.0               # Lattice spacing
#    self.K      = 1.165             # Coefficient ot solve Fermion doubling
#    self.gamma  = 3.16*1.42*3/2     # gamma = hbar * v_F
#    self.ptable = atoms.PeriodicTable()   # Periodic table for graphene k.p
#    self.ptable.add(0, "D", 2, 2)   # Fake atom for k.p
    
hamiltonian.GrapheneKpParams.__init__ = _GrapheneKpParamsInit

# Graphene k.p pickle support
def _GrapheneKpParamsSetState(self, dct):
    self.__dict__.update(dct)
    
hamiltonian.GrapheneKpParams.__setstate__ = _GrapheneKpParamsSetState

def _GrapheneKpParamsGetState(self):
    dct = dict(self.__dict__)
    return dct

hamiltonian.GrapheneKpParams.__getstate__ = _GrapheneKpParamsGetState
hamiltonian.GrapheneKpParams.__getstate_manages_dict__ = True


# Graphene tight binding default parameters
_GrapheneTbParamsOrgInit =  hamiltonian.GrapheneTbParams.__init__
def _GrapheneTbParamsInit(self):
    _GrapheneTbParamsOrgInit(self)
#    self.dtol   = 1E-3              # Tolerance when considering neighbors
#    self.ec         = 0             # Site energy of C
#    self.di0        = 1.42          # In-plane C-C bond length 
#    self.ti0        = 3.16          # In plane C-C hopping parameter
#    self.do0        = 3.35          # Out-of-plane C-C distance
#    self.to0        = 0.39          # Out-of-pane C-C hopping parameter
#    self.lmdz       = 0.6           # For interlayer all nearest neighbors
#    self.lmdxy      = 1.7           # See PRL 109, 236604 (2012)
#    self.alpha      = 1.65
#    self.doX        = 6.0           # Out-of-plane neighbor cut-off distance
#    self.ptable = atoms.PeriodicTable()   # Periodic table for graphene
    
hamiltonian.GrapheneTbParams.__init__ = _GrapheneTbParamsInit

# Graphene tight binding pickle support
def _GrapheneTbParamsSetState(self, dct):
    self.__dict__.update(dct)
    
hamiltonian.GrapheneTbParams.__setstate__ = _GrapheneTbParamsSetState

def _GrapheneTbParamsGetState(self):
    dct = dict(self.__dict__)
    return dct

hamiltonian.GrapheneTbParams.__getstate__ = _GrapheneTbParamsGetState
hamiltonian.GrapheneTbParams.__getstate_manages_dict__ = True

# Coherent RGFA default parameters
#_CohRgfaParamsOrgInit =  negf.CohRgfaParams.__init__
#def _CohRgfaParamsInit(self, nb):
#    _CohRgfaParamsOrgInit(self, nb)
#    self.kT             = 0.0259        # Temperature in eV (300 K)
#    self.ieta           = 1E-3j         # Contact imaginary potential
#    self.mu             = 0.0           # Device Fermi level
#    self.muS            = 0.0           # Source Fermi level
#    self.muD            = 0.0           # Source Fermi level
#    self.isOrthogonal   = True          # Orthogonal basis?
#    self.Emin           =-1.0           # Minimum energy 
#    self.Emax           = 1.0           # Maximum energy
#    self.dE             = 0.005         # Energy step
#    self.AutoGenE       = False         # Generate grid automatically?
#    self.AdaptiveGrid   = False         # Adaptive E grid?
    
# negf.CohRgfaParams.__init__ = _CohRgfaParamsInit

#_CohRgfaParamsOrgStr =  negf.CohRgfaParams.__str__
#def _CohRgfaParamsStr(self):
#    msg = _CohRgfaParamsOrgStr(self)
#    msg += "\n" + self.Prefix + " Emin         = " + str(self.Emin)
#    msg += "\n" + self.Prefix + " Emax         = " + str(self.Emax)
#    msg += "\n" + self.Prefix + " dE           = " + str(self.dE)
#    msg += "\n" + self.Prefix + " AutoGenE     = " + str(self.AutoGenE)
#    msg += "\n" + self.Prefix + " AdaptiveGrid = " + str(self.AdaptiveGrid)
#    
#    return msg
#negf.CohRgfaParams.__str__ = _CohRgfaParamsStr

# Pickle support
#def _CohRgfaParamsSetState(self, dct):
#    self.__dict__.update(dct)
#    
#negf.CohRgfaParams.__setstate__ = _CohRgfaParamsSetState
#
#def _CohRgfaParamsGetState(self):
#    dct = dict(self.__dict__)
#    return dct
#
#negf.CohRgfaParams.__getstate__ = _CohRgfaParamsGetState
#negf.CohRgfaParams.__getstate_manages_dict__ = True


#def _CohRgfaParamsGetInitArgs(self):
#    return (self.nb,)

#negf.CohRgfaParams.__getinitargs__ = _CohRgfaParamsGetInitArgs

# Electrostatic potential
_LinearPotOrgInit =  potential.LinearPot.__init__
def _LinearPotInit(self, atoms = None):
    _LinearPotOrgInit(self, atoms)
    self.rVS        =-0.5          # Fraction of voltage to be applied to drain
    self.rVD        = 0.5          # Fraction of voltage to be applied to drain
    self.rVG        = []           # Gate voltage ratios for VGG
    self.rVo        = []           # Gate voltage ratios for Vo
    self.gql        = []           # Gate quadrilaterals just for pickle
    self.lql        = []           # Linear region quadrilaterals just for pickle

potential.LinearPot.__init__ = _LinearPotInit

# Pickle support
def _LinearPotSetState(self, dct):
    self.__dict__.update(dct)
    
potential.LinearPot.__setstate__ = _LinearPotSetState

def _LinearPotGetState(self):
    dct = dict(self.__dict__)
    return dct

potential.LinearPot.__getstate__ = _LinearPotGetState
potential.LinearPot.__getstate_manages_dict__ = True

##
# Enum pickler.
#
def isEnumType(o):
    return isinstance(o, type) and issubclass(o,int) and not (o is int)

def _tuple2enum(enum, value):
    enum = getattr(utils, enum)
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
    for e in [ e for e in vars(utils).itervalues() if isEnumType(e) ]:
        pickle(e, reduce_enum)

_registerEnumPicklers()



