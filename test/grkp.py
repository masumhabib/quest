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
        self.workers = workers
        self.masterID = 0
        self.myID = self.workers.rank
        self.iAmMaster = (self.myID == self.masterID)
        
        # simulation 
        self.VDmin      = 0
        self.VDmax      = 0
        self.dVD        = 0.05
        self.VGmin      = 0
        self.VGmax      = 0
        self.dVG        = 0.05

        # output settings
        self.OutFileName = "gnr"
        
        # Device parameters
        self.nl         = 5
        self.nw         = 3

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
        
        
    def run(self):
        pass
 

def simulate(workers):

    d = GrapheneTransport(workers)
    d.prepare()

    """# prepare device
    p = GrapheneKpParams()
    p.ds    = 1
    p.rVD   = 0.5
    p.rVS   = -0.5
    p.rVG0  = 1.0
    p.rVG1  = 1.0
    #d = Device( p)
    d = Device(workers, p)

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
