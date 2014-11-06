#!/usr/bin/python

"""
 TI band structure simulator using QMICAD.
 
 Author: K M Masum Habib <masum.habib@virginia.edu>
 Last update: 05/08/2014
"""

import  sys
import  numpy as np
from    math import pi
import  mpi

import  qmicad
from    qmicad.simulators.dirackp import * 


##
# Run the simulation
def simulate(workers):

    # some constants 
    # Transport simulator
    bs = Band(workers)
    bs.HamType = bs.HAM_TI_SURF_KP
    bs.verbosity = vprint.MSG_NORMAL
    
    #Simulation parameters.
    bs.VERSION_REQUIRED = "0.08.2"

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
        str = "QMICAD version mismatch. Required "
        str += str(bs.VERSION_REQUIRED) + ". Got "
        str += str(qmicad.version)
        raise(str)

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

    # get MPI communicator
    workers = mpi.world
    # Run the simulator
    simulate(workers)
        
    return 0

"""
 Entry point.
"""
if __name__ == "__main__":
    sys.exit(main())

