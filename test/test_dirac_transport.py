#!/usr/bin/python

"""
 TI Surface transport simulator for p-n junction using QMICAD.
 
 Author: K M Masum Habib <masum.habib@virginia.edu>
 Last update: 03/26/2014
"""

import  sys
import  numpy as np
import  math
from    math import pi, tan, cos, sin

import  mpi
import  qmicad
from    qmicad.simulators.utils.linspace import linspace
from    qmicad.simulators.dirackp import *


##
# Run the simulation
def simulate(workers):

    # some constants 
    # Transport simulator
    tr = Transport(workers)
    tr.HamType = tr.HAM_TI_SURF_KP
    tr.verbosity = vprint.MSG_NORMAL
    
    # Simulation parameters ------------------------------------------
    tr.VERSION_REQUIRED = "0.06.5"
    tr.verbosity = vprint.MSG_DEBUG 
    # Do a dry run    
    tr.DryRun = False
    #tr.DryRun = True

    # Device structure  
    tr.gates       = [];         # Gate tuple (width, theta, voltage ratio)
    tr.a           = 5
    tr.K           = 1

    # Device structure
    tr.nw         = 21 
    tr.th         = 0
    tr.ng1        = 10 
    tr.ns1        = 0          
    tr.ng2        = 10
    
    thr = tr.th*pi/180;
    # gate # 1 right edge
    tr.addGateEdge(tr.ng1, 0, 1.0, 0.0)            
    # split # 1 right edge
    tr.addGateEdge(tr.ns1, 0)            
    # gate # 2 right edge
    tr.addGateEdge(tr.ng2 + tr.nw/2*tan(thr), tr.th, -1.0, 0.0)            
 
    # Get the total length
    tr.nb = tr.computeDevLen();

    # Output path
    tr.OutPath     = "./pn_TE"

    # Calculation type
    tr.Calculations["TE"] = 1
    #tr.Calculations["n"] = [{"N":2, "Block":'All'}]
    tr.Calculations["n"] = [{"N":2, "Block":4}]
    #tr.Calculations["I"] = [{"N":2, "Block":'All'}]
    tr.Calculations["I"] = [{"N":2, "From":0, "To":1},
                            {"N":2, "From":tr.nb-2, "To":tr.nb-1}]


    # Output path
    tr.OutPath    += "/"
    tr.OutFileName = "TR"
    

    # Bias
    tr.VGG        = np.array([0.0])    # Gate voltage offset
    tr.Vo         = 0.15               # Built-in voltage

    tr.VDD        = np.array([0.0])
    tr.rVS        =-0.5       # source ratio
    tr.rVD        = 0.5       # drain ratio

    # Energy range
    tr.Emin       =-0.5             # Minimum energy 
    tr.dE         = 0.005            # Energy step
    tr.Emax       = 0.5-tr.dE+0.00001       # Maximum energy

    tr.kT         = 0.0259  # temperature

    # --------------------------------------------------------------
    if qmicad.version != tr.VERSION_REQUIRED:
        msg = "QMICAD version mismatch. Required "
        msg += str(tr.VERSION_REQUIRED) + ". Got "
        msg += str(qmicad.version)
        raise Exception(msg)

    
    # Get the total length
    tr.nb = tr.computeDevLen()

    # Change lattice constant if requested
    if hasattr(tr, 'a'):
        if not hasattr(tr, 'hp'):
            tr.hp = TISurfKpParams()
        tr.hp.a = tr.a
    
    if hasattr(tr, 'K'):
        if not hasattr(tr, 'hp'):
            tr.hp = TISurfKpParams()
        tr.hp.K = tr.K
 
    # Create the atomistic geometry
    tr.createAtomicGeom()
    
    # Add roughness if requested
    if hasattr(tr, 'roughness'):    
        tr.createRoughEdges(tr.roughness)
    

    # Generate the Hamiltonian
    tr.generateHamiltonian()
    
    # Set up the potential solver
    tr.setupPotential()
    
    #
    # Gates and contacts
    #
    xmn = tr.xmn
    xmx = tr.xmx
    ymn = tr.ymn
    ymx = tr.ymx
    a  = tr.hp.a
    
    # Source contact    
    beg = xmn
    end = beg + a
    tmp = Quadrilateral(Point(beg, ymn), Point(end, ymn),
                        Point(end, ymx), Point(beg, ymx))

    if hasattr(tr, 'rVS'):
        tr.addSource(tmp, tr.rVS)
    else:
        tr.addSource(tmp)

    # Gates and splits
    beg = end
    begb = end
    begt = end
    for gate in tr.gates:
        tanth = tan(gate["th"]*pi/180.0)
        sd    = tr.nw*a/2.0*tanth
        end   = beg + gate["nl"]*a
        endb  = end - sd
        endt  = end + sd
        
        tmp = Quadrilateral(Point(begb, ymn),  Point(endb, ymn),
                            Point(endt, ymx),  Point(begt, ymx))
        
        beg = end
        begb = endb
        begt = endt

        rVo = gate["rVo"]
        rVG = gate["rVG"]
        if rVo is None and rVo is None:
            tr.addLinearRegion(tmp)
        else:
            tr.addGate(tmp, rVo, rVG)
            
    # Drain contact
    begb = endb
    begt = endt
    endb = begb + a
    endt = begt + a
    tmp = Quadrilateral(Point(begb, ymn), Point(xmx, ymn),
                        Point(xmx, ymx), Point(begt, ymx))
    if hasattr(tr, 'rVD'):
        tr.addDrain(tmp, tr.rVD)
    else:
        tr.addDrain(tmp)

            
    # Run
    tr.run()

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

# Inject some helpers function to Transport
def _addGateEdge(self, nl, th, rVo=None, rVG=None):
    """Adds gate or a linear region."""
    self.gates.append({"nl":int(nl), "th":th, "rVo":rVo, "rVG":rVG})
Transport.addGateEdge = _addGateEdge

def _computeDevLen(self):
    nb = 3;
    for gate in self.gates:
        nb = nb + gate["nl"]
    return int(nb)
Transport.computeDevLen = _computeDevLen


"""
 Entry point.
"""
if __name__ == "__main__":
    sys.exit(main())



