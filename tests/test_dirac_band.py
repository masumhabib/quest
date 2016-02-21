#!/usr/bin/python

"""
 TI band structure simulator using QMICAD.
 
 Author: K M Masum Habib <masum.habib@virginia.edu>
 Last update: 05/08/2014
"""

import  sys
import  numpy as np
from    math import pi

# Workaround for a bug involving Boost.MPI, OpenMPI and Python
# in linux system. You can ignore these lines if you are not using
# OpenMPI in linux.
#import sys
#if sys.platform == 'linux2':
#    import DLFCN as dl
#    flags = sys.getdlopenflags()
#    sys.setdlopenflags(dl.RTLD_NOW|dl.RTLD_GLOBAL)
#    import mpi
#    sys.setdlopenflags(flags)
#else:
#    import mpi

import  qmicad
from    qmicad.simulators.dirackp import * 


##
# Run the simulation
def simulate():

    # some constants 
    # Band structure simulator
    bs = Band()
    bs.HamType = bs.HAM_TI_SURF_KP
    bs.verbosity = vprint.MSG_NORMAL
    
    #Simulation parameters.
    bs.VERSION_REQUIRED = "0.10.1"

    bs.DryRun = False
    bs.AutoGenKp = False

    # Cell structure
    bs.nw         = 1         # width
    bs.nl         = 1         # length
    bs.nh         = 1         # height

    bs.dk         = 0.05       
    bs.Dim        = 2         # 2D structure.
    bs.a          = 5

    bs.nk1        = 101
    bs.nk2        = 101 
    bs.nb         = 4                # Number of bands to save

    # Output path
    bs.OutPath     = "./2D"

    # Output path
    bs.OutPath    += "/"
    bs.OutFileName = "EK"

    lx = bs.nl*bs.a;
    bs.kp.addKLine(Point(-pi/lx, 0), Point(pi/lx, 0), bs.nk1)

    # Setup and simulation
    if qmicad.version != bs.VERSION_REQUIRED:
        msg = "QMICAD version mismatch. Required "
        msg += str(bs.VERSION_REQUIRED) + ". Got "
        msg += str(qmicad.version)
        raise Exception(msg)

    # Change lattice constant if requested
    if hasattr(bs, 'a'):
        if bs.HamType == bs.HAM_TI_SURF_KP:
            bs.hp = TISurfKpParams()
            bs.hp.a = bs.a
    
    # Create the atomistic geometry
    bs.createAtomicGeom()
    
    # Create neighboring cells
    bs.addNeighbor()
    
    # Generate the Hamiltonian
    bs.generateHamiltonian()
        
    # Automatic generation of k points
    if bs.AutoGenKp == True:
        bs.generateKPoints()
        
                
    # Run
    bs.run()

"""
 The main() function.
"""
def main(argv = None):
    if argv is None:
        argv = sys.argv

    # Run the simulator
    simulate()
        
    return 0

"""
 Entry point.
"""
if __name__ == "__main__":
    sys.exit(main())

