"""
 Quantum Mechanics Inspired Computer Aided Design (QMICAD) library is a 
 collection of C++ classes and functions for simulation and design of 
 nano-scaled devices using quantum mechanical tools such as Non-equilibrium 
 Green's Function (NEGF) formalism. This library provides a common framework 
 for NEGF so that it can be used with any empirical tight binding, k.p model, 
 extended Huckel method and density functional theory codes. This library also
 includes a generic empirical tight binding model and a k.p model which can be 
 extended for any material with known parameters.
  
 Author: K M Masum Habib <masum.habib@virginia.edu>
 Last update: 02/12/2014

"""
from _qmicad import version 
from _qmicad import greet
from _qmicad import setVerbosity

import utils
import atoms
import kpoints
import potential
import hamiltonian
import band
import negf
import simulators