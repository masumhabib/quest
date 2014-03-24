""" 
    Non-equilibrium Green Function Formalism simulation package. It contains the 
    following:
        * Coherent Recursive Green Function Algorithm (RGFA)
        * Parallel Energy Loop
    
    Author: K M Masum Habib <masum.habib@virginia.edu>
    Last update: 02/12/2014
"""
 
from _negf import CohRgfaParams
from _negf import NegfEloop

# Coherent RGFA default parameters
_CohRgfaParamsOrgInit =  CohRgfaParamsParams.__init__
def _CohRgfaParamsParamsInit(self):
    _CohRgfaParamsOrgInit(self)
    self.kT             = 0.0259        # Temperature in eV (300 K)
    self.eta            = 1E-3          # Contact imaginary potential
    self.mu             = 0.0           # Device Fermi level
    self.isOrthogonal   = True          # Orthogonal basis?
    self.Emin           =-1.0           # Minimum energy 
    self.Emax           = 1.0           # Maximum energy
    self.dE             = 0.005         # Energy step
    self.AutoGenE       = False         # Generate grid automatically?
    self.AdaptiveGrid   = False         # Adaptive E grid?
    