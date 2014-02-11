# Graphene k.p simulator using QMICAD
# Author: K M Masum Habib <masum.habib@virginia.edu>
#
#
#
from ctypes import *
import sys
import getopt
import mpi
from qmicad import *

class GrapheneTransport:
    def __init__(self, workers):        
        # MPI stuff
        self.workers = workers                          # MPI workers
        self.masterID = 0                               # Master ID
        self.myID = self.workers.rank                   # This process ID
        self.iAmMaster = (self.myID == self.masterID)   # Master/Slave?
        
        # Bias
        self.VDmin      = 0             # Max drain-source bias
        self.VDmax      = 0             # Min drain-source bias
        self.dVD        = 0.05          # Drain-source bias step
        self.VGmin      = 0             # Min Gate bias
        self.VGmax      = 0             # Max Gate bias
        self.dVG        = 0.05          # Gate bias step
        self.rVD        = 0.5           # Fraction of voltage to be applied to drain
        self.rVS        = -0.5          # Fraction of voltage to be applied to drain
        self.rVG1       = 1.0           # Fraction of voltage to be applied to gate # 1
        self.rVG2       = 1.0           # Fraction of voltage to be applied to gate # 2

        # output settings
        self.OutFileName = "gnr"        # Output file prefix
        
        # Device parameters
        self.nl         = 5             # Length of the device+contacts
        self.nw         = 3             # Width of the device
        self.ds         = 1             # Split width: distance between gate 1 and 2.

        # NEGF parameters
        self.kT         = 0.0259
        self.eta        = 1E-3
        self.mu         = 0.0
        self.Emin       = -3.1
        self.Emax       = 3.1
        self.dE         = 0.1
        self.AutoGenE   = False

        # graphene k.p parameters
        self.ax         = 4
        self.ay         = 4
        self.dtol       = 1E-3
        self.K          = 1
        self.gamma      = 3.16*1.42*3/2

        # Periodic table with fake atoms
        self.ptable     = PeriodicTable()
        self.ptable.add(0, "D", 2, 2)

  
    ## 
    # Prepare the sumulation.
    def prepare(self):
        # Print welcome message.
        if( self.iAmMaster):
            print greet()
            
        # create atomistic geometry of the device.
        self.geom = AtomicStruct(self.nl, self.nw, self.ax, self.ay, self.ptable)
        
        # k.p parameters for graphene
        self.grkpp = GrapheneKpParams()
        self.grkpp.dtol = self.dtol
        self.grkpp.ax = self.ax
        self.grkpp.ay = self.ay
        self.grkpp.K = self.K
        self.grkpp.gamma = self.gamma
        self.grkpp.update()
        
        # generate hamiltonian and overlap matrices
        self.ham = GrapheneKpHam(self.grkpp)            # hamiltonian generator
        lyr0 = self.geom.span(0, self.nw-1);            # extract block # 0
        lyr1 = self.geom.span(self.nw, 2*self.nw-1);    # extract block # 1
        self.ham.setSize(2)                             # uniform device, only need H_0,0 and H_1,0
        self.ham.genDiagBlock(lyr0, lyr0, 0)            # generate H_0,0 and S_0,0
        self.ham.genLowDiagBlock(lyr1, lyr0, 0)         # generate H_1,0 
        
        # Setup the NEGF parameters
        self.np = NegfParams(self.nl)
        for ib in range(0, self.nl+1):                  # setup the block hamiltonian
            if (ib != self.nl):
                self.np.H0(self.ham.H0(0), ib)          # H0: 0 to N+1
                self.np.S0(self.ham.S0(0), ib)          # S0: 0 to N+1
            self.np.Hl(self.ham.Hl(0), ib)              # Hl: 0 to N+2
        self.np.isOrthogonal = True                     # k.p is an orthogonal basis
        self.np.kT = self.kT;                           
        self.np.ieta = 1j*self.eta;                 # imaginary potential
        self.np.grcCache = Option.Enabled;          # enable grc cache
        self.np.DCache = Option.Enabled;            # enable D cache

        # Electrostatic potential
    """# prepare device
    p = GrapheneKpParams()

    if (iAmMaster):
        d.tic()

    d.prepare()

    # add contacts
    xmn = d.xmin() - p.ax/2
    xmx = d.xmax() + p.ax/2
    ymn = d.ymin() - p.ay/2
    ymx = d.ymax() + p.ay/2
    ds = (p.ds-1)*p.ax;

    # source
    ql = Quadrilateral(Point(xmn, ymn), Point(xmn+p.ax, ymn),
                       Point(xmn+p.ax, ymx),  Point(xmn, ymx))
    d.addSource(ql)
    # drain
    ql = Quadrilateral(Point(xmx-p.ax, ymn), Point(xmx, ymn),
                       Point(xmx, ymx), Point(xmx-p.ax, ymx))
    d.addDrain(ql)
    # gates
    ql = Quadrilateral(Point(xmn+p.ax, ymn), Point(-ds/2, ymn), 
                       Point(-ds/2, ymx), Point(xmn+p.ax, ymx))
    d.addGate(ql)
    ql = Quadrilateral(Point(ds/2, ymn), Point(xmx-p.ax, ymn),
                       Point(xmx-p.ax, ymx), Point(ds/2, ymx))
    d.addGate(ql)
    # linear region
    ql = Quadrilateral(Point(-ds/2, ymn), Point(ds/2, ymn), 
                       Point(ds/2, ymx), Point(-ds/2, ymx))
    d.addLinearRegion(ql)
"""        
        
    def run(self):
        pass
 

def simulate(workers):

    d = GrapheneTransport(workers)
    d.prepare()


"""
    # bias loops
    # loop over drain
    for ivd in range(0, d.NVDD()):
        VDD = d.VDD(ivd);
        d.VDS(VDD*p.rVD, VDD*p.rVS)
        # loop over gate
        for ivg in range(0, d.NVGG()):
            VGG = d.VGG(ivg)
            d.VG(0, VGG*p.rVG0)
            d.VG(0, VGG*p.rVG0)
            d.VLR(0, VGG*p.rVG0, VGG*p.rVG1)
            d.computePotential()
            np = d.NegfParam()
            Eloop = QneffEloop(np, workers)
            d.runNegfEloop()

    #print d.toString()
    
    
    if(iAmMaster):
        d.toc()
        print d.time()
"""

# main function
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # get MPI communicator
    workers = mpi.world
    # Run the simulator
    simulate(workers)
    return 0


if __name__ == "__main__":
    sys.exit(main())
