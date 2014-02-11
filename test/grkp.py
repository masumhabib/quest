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
        grkpp = GrapheneKpParams()
        grkpp.dtol = self.dtol
        grkpp.ax = self.ax
        grkpp.ay = self.ay
        grkpp.K = self.K
        grkpp.gamma = self.gamma
        grkpp.update()
        
        # generate hamiltonian and overlap matrices
        self.ham = GrapheneKpHam(grkpp)                 # hamiltonian generator
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
        self.np.ieta = 1j*self.eta;                     # imaginary potential
        self.np.grcCache = Option.Enabled;              # enable grc cache
        self.np.DCache = Option.Enabled;                # enable D cache

        # Electrostatic potential
        self.V = LinearPot(self.geom)
        # add contacts
        xmn = self.geom.xmin() - self.ax/2
        xmx = self.geom.xmax() + self.ax/2
        ymn = self.geom.ymin() - self.ay/2
        ymx = self.geom.ymax() + self.ay/2
        ds = (self.ds-1)*self.ax;
        # source
        ql = Quadrilateral(Point(xmn, ymn), Point(xmn+self.ax, ymn),
                           Point(xmn+self.ax, ymx),  Point(xmn, ymx))
        self.V.addSource(ql)
        # drain
        ql = Quadrilateral(Point(xmx-self.ax, ymn), Point(xmx, ymn),
                           Point(xmx, ymx), Point(xmx-self.ax, ymx))
        self.V.addDrain(ql)
        # gates
        ql = Quadrilateral(Point(xmn+self.ax, ymn), Point(-ds/2, ymn), 
                           Point(-ds/2, ymx), Point(xmn+self.ax, ymx))
        self.V.addGate(ql)
        ql = Quadrilateral(Point(ds/2, ymn), Point(xmx-self.ax, ymn),
                           Point(xmx-self.ax, ymx), Point(ds/2, ymx))
        self.V.addGate(ql)
        # linear region
        ql = Quadrilateral(Point(-ds/2, ymn), Point(ds/2, ymn), 
                           Point(ds/2, ymx), Point(-ds/2, ymx))
        self.V.addLinearRegion(ql)
        
        #Bias grid
        self.vdd = VecGrid(self.VDmin, self.VDmax, self.dVD)
        self.vgg = VecGrid(self.VGmin, self.VGmax, self.dVG)
        
        
    def run(self):
        # Bias loops
        # loop over drain
        for ivd in range(0, self.vdd.N()):
            VDD = self.vdd.V(ivd)
            VD = VDD*self.rVD
            VS = VDD*self.rVS
            self.np.muS = self.mu - VS
            self.np.muD = self.mu - VD
            # loop over gate
            for ivg in range(0, self.vgg.N()):
                VGG = self.vgg.V(ivg)
                VG1 = VGG*self.rVG1
                VG2 = VGG*self.rVG1
                self.V.VG(0, VG1)               # bias voltage for gate#1 
                self.V.VG(1, VG1)               # bias voltage for gate#2
                self.V.VS(VG1)                  # potential of source
                self.V.VD(VG2)                  # potential of drain
                self.V.VLR(0, VG1, VG2)         # potential of linear region
                # solve for potential
                self.V.compute()
                
                # export potential to NEGF
                for ib in range(0, self.nl):                  # setup the block hamiltonian
                    self.np.V(self.V.toOrbPot(self.nw*ib, self.nw*(ib+1)-1), ib)

                # create energy grid
                if (self.AutoGenE):
                    Emin = self.np.muD - 10*self.kT
                    Emax = self.np.muS + 10*self.kT
                else:
                    Emin = self.Emin
                    Emax = self.Emax
                EE = VecGrid(Emin, Emax, self.dE)
                
                # create and run energy loop
                Eloop = NegfEloop(EE, self.np, self.workers) 
                Eloop.run()
                
                # save results
                if (self.iAmMaster):
                    fileName = self.OutFileName + "TE_VDD" + str(VDD) + "VGG" + str(VGG) + ".dat"
                    Eloop.saveTE(fileName)
##
# Run the simulation
def simulate(workers):
    # Start clock
    clock = Timer()
    clock.tic()

    gt = GrapheneTransport(workers)
    gt.prepare()
    gt.run()
    
    clock.toc()
    print clock

##
# main function.
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # get MPI communicator
    workers = mpi.world
    # Run the simulator
    simulate(workers)
    return 0

##
# Entry point.
if __name__ == "__main__":
    sys.exit(main())
