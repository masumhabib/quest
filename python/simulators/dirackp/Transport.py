"""
    Qunantum Transport simulator using QMICAD.
    Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>

    Last update: 05/12/2014
"""

import os
import pickle as pk
import numpy as np
import random as rn

from qmicad import setVerbosity, greet
from qmicad.atoms import AtomicStruct, SVec, LCoord
from qmicad.hamiltonian import TISurfKpParams4, TISurfKpParams, TI3DKpParams, GrapheneKpParams, GrapheneTbParams, generateHamOvl
from qmicad.negf import CohRgfLoop
from qmicad.kpoints import KPoints
from qmicad.potential import LinearPot
from qmicad.utils import Timer, Workers, Quadrilateral, Point
from qmicad.simulators.utils import vprint
from qmicad.simulators.utils.vprint import nprint, dprint, eprint
from qmicad.simulators.utils.linspace import linspace

class Transport(object):
    """
    Qunantum Transport simulator using QMICAD.
        
    * Device geometry:
        ----------------------------------------
         ... |-1 | 0 | 1 |  ...  | N |N+1|N+2| ...
        ----------------------------------------
                   ^  <------^------>  ^
                 left      Device    right
               contact              contact

    Notes:
        1. Coherent RGF calculation for non-uniform devices requires at least 
           5 blocks: 1 source block, 1 drain block and 3 device blocks as 
           shown in the following fig:

           --------------------------
             | 0 | 1 | 2 | 3 | 4 |
           --------------------------
               ^  <--------->  ^
            source  Device   drain 

           Block  0 and block 1 has to be identical. 
           Simularly, block 3 and block 4 hast to be identical.

        2. Any uniform device needs atleast 3 blocks: two contact
           blocks and one device block.
           
    Attributes:
    
        
    """

    def __init__(self, workers):   
        """ 
        Class constructor. 
        param:
          workers: mpi.world - MPI world communicator.
        """  

        self.verbosity = vprint.MSG_NORMAL # Verbosity level
        
        # MPI stuff
        self.workers = Workers(workers)
        vprint.IAmMaster = self.workers.IAmMaster()
    
        # Timer
        self.clock = Timer()

        # output path settings
        self.OutPath    = "./out/"     # Output path
        self.OutFileName = "TR"        # Output file name prefix
        
        # Bias
        self.VDD        = np.zeros(1)  # Drain bias
        self.VGG        = np.zeros(1)  # Gate bias
        self.Vo         = 0.0          # Built in potential
        
        # Device geometry
        self.nb         = 11           # Length of the device+contacts
        self.nw         = 9            # Width of the device
        self.nh         = 1
        self.nbw        = []           # Width of all layers

        # for k-loop
        self.nk2         = 0            # number of k points along a2        
        self.kp         = KPoints()    # K point generator
        
        # Hamiltonian type
        self.HAM_TI_SURF_KP  = 10       # TI surface k.p hamiltonian
        self.HAM_TI_SURF_KP4 = 11       # TI surface k.p hamiltonian with 4 spin basis set
        self.HAM_TI_3D_KP    = 15       # TI 3D k.p hamiltonian with 4 spin basis set        
        self.HAM_GRAPHENE_KP = 20       # Graphene k.p Hamiltonian
        self.HamType         = self.HAM_TI_SURF_KP 
        
        # Device types
        self.COH_RGF_UNI    = 10       # Coherent uniform RGF type device
        self.COH_RGF_NON_UNI= 20       # Coherent Non-uniform RGF type device
        self.DevType        = self.COH_RGF_UNI
        
        # Potential profile type
        self.POT_LINEAR     = 10       # Linear potential model
        self.PotType        = self.POT_LINEAR  
        
        # Calculations
        self.Calculations = {}
        self.Calculations["TE"] = 1
        
        # 
        self.kT             = 0.0259        # Temperature in eV (300 K)
        self.ieta           = 1E-3j         # Contact imaginary potential
        self.mu             = 0.0           # Device Fermi level
        self.OrthoBasis     = True          # Orthogonal basis?
        self.Emin           =-1.0           # Minimum energy 
        self.Emax           = 1.0           # Maximum energy
        self.dE             = 0.005         # Energy step
        self.AutoGenE       = False         # Generate grid automatically?
        
        # Debug stuffs
        self.DebugPotFile   = "dbg_pot.dat"
        self.DebugGeomFile  = "dbg_geom.svg"
        
        # Dry run
        self.DryRun         = False
        
        # Skip if resulting dat file already exists?
        self.SkipExistingSimulation = False

        # Print welcome message.
        nprint(greet())

    @property 
    def verbosity(self):
        return self._verbosity  
    @verbosity.setter
    def verbosity(self, value):
        self._verbosity = value
        # Set verbosity level of QMICAD library.
        setVerbosity(self._verbosity)
        vprint.verbosity = self._verbosity
   
    def muS(self, VDD):
        """Calculates the fermi energy level of the source."""
        return self.mu - VDD*self.V.rVS
    
    def muD(self, VDD):
        """Calculates the drain energy level of the source."""
        return self.mu - VDD*self.V.rVD
 
    
    def VG(self, VGG, Vo, ig):
        """ Computes the gate voltage of gate #ig. """
        # Set gate voltages
        return VGG*self.V.rVG[ig] + Vo*self.V.rVo[ig]
      
    def addSource(self, sql, rVS = -0.5):
        """ Adds a source contact to the device for electrostatics calculation."""
        self.V.addSource(sql)
        self.V.rVS = rVS
        # just for pickling
        self.V.sql = sql    
    
    def addDrain(self, dql, rVD = 0.5):
        """ Adds a drain contact to the device for electrostatics calculation."""
        self.V.addDrain(dql)
        self.V.rVD = rVD
        # just for pickling
        self.V.dql = dql 

    def addGate(self, gql, rVo = 1, rVG = 1):
        """ Adds a gate to the device for electrostatics calculation."""
        self.V.addGate(gql)
        self.V.rVG.append(rVG)
        self.V.rVo.append(rVo)
         # just for pickling
        self.V.gql.append(gql) 

    def addLinearRegion(self, lql):
        """ Adds a linear region to the device for electrostatics calculation."""
        self.V.addLinearRegion(lql)
        # just for pickling
        self.V.lql.append(lql) 
               
    def createAtomicGeom(self):
        """ Creates atomistic geometry. """

        nprint("\n Creating atomistic geometry ...")

        # TI k.p surface      
        if (self.HamType == self.HAM_TI_SURF_KP):    
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = TISurfKpParams()
            else:
                if not isinstance(self.hp, TISurfKpParams): # if hp exists but not TISurfKpParams type, create it.
                    self.hp = TISurfKpParams()
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nb, self.nw, self.nh)
            self.nbw = [self.nw]*self.nb
        elif (self.HamType == self.HAM_TI_SURF_KP4):    
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = TISurfKpParams4()
            else:
                if not isinstance(self.hp, TISurfKpParams4): # if hp exists but not TISurfKpParams type, create it.
                    self.hp = TISurfKpParams4()
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nb, self.nw, self.nh)
            self.nbw = [self.nw]*self.nb
        elif (self.HamType == self.HAM_TI_3D_KP):    
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = TI3DKpParams()
            else:
                if not isinstance(self.hp, TI3DKpParams): # if hp exists but not TISurfKpParams type, create it.
                    self.hp = TI3DKpParams()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nb, self.nw, self.nh)
            self.nbw = [self.nw]*self.nb
        elif (self.HamType == self.HAM_GRAPHENE_KP):
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = GrapheneKpParams()
            else:
                if not isinstance(self.hp, GrapheneKpParams): # if hp exists but not GrapheneKpParams type, create it.
                    self.hp = GrapheneKpParams()
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nb, self.nw, self.nh)
            self.nbw = [self.nw]*self.nb
        else:
            raise RuntimeError(" Unsupported Hamiltonian type. ")
        
        self.updateBoundingBox()

        nprint(" done.")

    def createRoughEdges(self, sigma):
        """ 
        Creates rough edges. Works only for sorted lattice points.
        For rectangular lattice 
        """

        nprint("\n Creating rough edges ...")

        if (self.HamType == self.HAM_TI_SURF_KP 
                    or self.HamType == self.HAM_TI_SURF_KP4
                    or self.HamType == self.HAM_GRAPHENE_KP):
            nw = []
            geom = AtomicStruct()
            beg = 0
            # loop through the layers and 
            # remove some atoms randomly from the edges.
            for ib in range(0, self.nb):
                end = beg + self.nw - 1
                lyr = self.geom.span(beg, end)
                beg = end + 1
                
                if ( ib > 1 and ib < self.nb - 2):
                    # Remove or add some atoms from the bottom edge.
                    nr = int(rn.gauss(0, sigma))
                    while abs(nr) > 3*sigma:
                        nr = int(rn.gauss(0, sigma))
                    if nr <= 0: 
                        lyr = lyr.span(abs(nr), lyr.NumOfAtoms - 1)
                    else:
                        xtra = lyr.span(0, nr-1)
                        lv = SVec(0, - self.hp.ay*nr, 0)
                        xtra = xtra + lv
                        lyr = xtra + lyr 
                        
                    # Remove or add some atoms from the top edge
                    nr = int(rn.gauss(0, sigma))
                    while abs(nr) > 3*sigma:
                        nr = int(rn.gauss(0, sigma))
                    if nr <= 0: 
                        lyr = lyr.span(0, lyr.NumOfAtoms - 1 - abs(nr))
                    else:
                        xtra = lyr.span(lyr.NumOfAtoms - nr, lyr.NumOfAtoms - 1)
                        lv = SVec(0, self.hp.ay*nr, 0)
                        xtra = xtra + lv
                        lyr = lyr + xtra
                        
                # save this layer
                geom = geom + lyr
                nw.append(lyr.NumOfAtoms)
                
            # save geometry
            self.geom = geom
            self.nbw = nw
            self.DevType = self.COH_RGF_NON_UNI
            
            # update bounding box
            self.updateBoundingBox()

        else:           
            raise RuntimeError(" Unsupported Hamiltonian type. ")

        nprint(" done.")
                
    def generateHamiltonian(self):
        """ Generates hamiltonian and overlap matrices. """

        nprint("\n Generating Hamiltonian matrix ...")

        # For uniform RGF blocks
        if (self.DevType == self.COH_RGF_UNI): 
            # no k-loop, real space hamiltonian only
            if self.kp.N == 0:
                
                lyr0 = self.geom.span(0, self.nw*self.nh-1);               # extract block # 0
                lyr1 = self.geom.span(self.nw*self.nh, 2*self.nw*self.nh-1);    # extract block # 1                
                
                self.H0, self.S0 = generateHamOvl(self.hp, lyr0, lyr0)

                self.Hl, S = generateHamOvl(self.hp, lyr1, lyr0)

#                np.set_printoptions(linewidth=200)
#                print "\nS0\n"
#                print self.S0
#                print "\nH0\n"
#                print self.H0
#                print "\nHl\n"
#                print self.Hl
#                self.geom.exportGjf('dbg_geom.gjf')
#                lyr0.exportGjf('lyr0.gjf')
#                lyr1.exportGjf('lyr1.gjf')
#                lyr01 = lyr0+lyr1;
#                lyr01.exportGjf('lyr01.gjf')
                
            # nearest neighbors in transverse direction for k-loop
            else:
                self.H0 = [];
                self.S0 = [];
                self.Hl = [];
                self.pv = []; 
                self.pvl = []; 
                
                lv = self.geom.LatticeVector
                self.pv.append(lv*LCoord(0,0,0))
                self.pvl.append(lv*LCoord(0,0,0))
                lyr0 = self.geom.span(0, self.nw*self.nh-1);               # extract block # 0
                lyr1 = self.geom.span(self.nw*self.nh, 2*self.nw*self.nh-1);    # extract block # 1
                
                H, S = generateHamOvl(self.hp, lyr0, lyr0)
                self.H0.append(H); self.S0.append(S)
                H, S = generateHamOvl(self.hp, lyr1, lyr0)
                self.Hl.append(H)
                
                self.pv.append(lv*LCoord(0,1,0))
                lyr0top = lyr0 + self.pv[1]                # top neighbor of layer 0
                H, S = generateHamOvl(self.hp, lyr0, lyr0top)
                self.H0.append(H)
                
                self.pv.append(lv*LCoord(0, -1, 0))
                lyr0bot = lyr0 + self.pv[2]                # bottom neighbor of layer 0
                H, S = generateHamOvl(self.hp, lyr0, lyr0bot)
                self.H0.append(H)
               
                self.pvl.append(lv*LCoord(-1,1,0))
                H, S = generateHamOvl(self.hp, lyr1, lyr0top)
                self.Hl.append(H)
                
                self.pvl.append(lv*LCoord(-1,-1,0))
                H, S = generateHamOvl(self.hp, lyr1, lyr0bot)        
                self.Hl.append(H)
        # Non-uniform RGF blocks        
        elif (self.DevType == self.COH_RGF_NON_UNI):
            self.H0 = [];
            self.Hl = [];
            beg = 0
            for ib in range(0, self.nb):                    # setup the block hamiltonian
                end = beg + self.nbw[ib] - 1
                
                # generate H_i,i and S_i,i
                lyri = self.geom.span(beg, end);            # extract block # i             
                H0,S0 = generateHamOvl(self.hp, lyri, lyri)
                self.H0.append(H0)

                # generate H_i,i-1
                if (ib > 0):
                    #self.ham.genLowDiagBlock(lyri, lyrim1, ib) 
                    Hl,Sl = generateHamOvl(self.hp, lyri, lyrim1)
                    self.Hl.append(Hl)
                    #self.np.Hl(self.Hl(ib), ib)              # Hl: 1 to N+1=nb-1                
                else:
                    self.S0 = S0
                
                lyrim1 = lyri
                beg = end + 1
            
        nprint(" done.")
        
    def setupPotential(self):
        """ Sets up the potential profile """

        nprint("\n Setting up potential ...")
        # Linear
        if (self.PotType == self.POT_LINEAR):
            self.V = LinearPot(self.geom)
        nprint(" done.")
    
    def check(self):
        ret = True
        if (self.DevType == self.COH_RGF_UNI):
            ret = ret and (self.nb >= 3)
        if (self.DevType == self.COH_RGF_NON_UNI):
            ret = ret and (self.nb >= 5)
            
        return ret

    def runBiasStep(self, VGG, Vo, VDD):            
        """Runs the sumulation."""

        # Set drain and Fermi levels
        self.rgf.mu(self.muD(VDD), self.muS(VDD))
            
        nprint("\n Bias loop:")                                
        nprint("\n  VGG = " + str(VGG) + ", Vo = " + str(Vo) 
                    + ", VDD = " + str(VDD) + ".")
        
        fileName = self.OutFileName + "_VGG{0:2.3f}_Vo{1:2.3f}_VDD{2:2.3f}".format(VGG, Vo, VDD)

        # skip calculation if result file exists.
        if self.SkipExistingSimulation == True:
            if os.path.isfile(self.OutPath + fileName + ".dat"):
                nprint("\n  Result exists, skipping.")
                return
            
        # Set gate voltages
        for ig in range(self.V.NG):
            VG = self.VG(VGG, Vo, ig)
            # bias voltage for gate # ig
            self.V.VG(ig, VG)
            # For linear potential profile
            if (self.PotType == self.POT_LINEAR):
                # set potential of source
                if (ig == 0):
                    self.V.VS(VG)
                # set potential of drain
                if (ig == self.V.NG - 1):
                    self.V.VD(VG)
        
        # Potential of linear region
        if (self.PotType == self.POT_LINEAR):
            for il in range(self.V.NLR):
                VG1 = self.VG(VGG, Vo, il)
                VG2 = self.VG(VGG, Vo, il+1)
                self.V.VLR(il, VG1, VG2)        

        # Print some useful information
        nprint("\n")
        nprint(self.dynamicnstr())
 
        # Compute the electrostatic potential
        self.V.compute()
        eprint(self.V)
        
        # Save the calculated potential.
        if (self.verbosity >= vprint.MSG_DEBUG):
            if (self.workers.IAmMaster()):
                # Create directory if not exist.
                if not os.path.exists(self.OutPath):
                    os.makedirs(self.OutPath) 
                self.V.exportPotential(self.OutPath + self.DebugPotFile)
        
        # Export potential to NEGF
        beg = 0
        self.Vo = []
        for ib in range(self.nb):                  # setup the block hamiltonian
            end = beg + self.nbw[ib]*self.nh - 1
            self.Vo.append(self.V.toOrbPot(beg, end)) 
            self.rgf.V(self.Vo[ib], ib)
            beg = end + 1
            
            
        # Create energy grid
        if (self.AutoGenE):
            Emin = self.muD(VDD) - 10*self.kT
            Emax = self.muS(VDD) + 10*self.kT
        else:
            Emin = self.Emin
            Emax = self.Emax
            
        EE = linspace(Emin, Emax, self.dE)

        
        nprint("\n Energy grid:\n  Min: " + str(Emin) + ", max: " 
                    + str(Emax) + ", interval " + str(self.dE)
                    + ".")
        nprint("\n Total " + str(len(EE)) + " energy point(s) " 
                    + "running on " + str(self.workers.N()) + " CPU(s): " 
                    + str(int(round(len(EE)/self.workers.N()))) + " pts/CPU ... \n")
                                
        # Set energy
        self.rgf.E(EE)

        # Run the simulation
        if (self.DryRun == False):
            self.rgf.run()
        nprint("\n  done.")
        
        nprint("\n Saving results to disk ...")
        # save results
        if (self.workers.IAmMaster()):
            # Create directory if not exist.
            if not os.path.exists(self.OutPath):
                os.makedirs(self.OutPath)            
            # save results to file.
            if (self.DryRun == False):
                fo = open(self.OutPath + fileName + ".dat", "wt")
                fo.close()
                self.rgf.save(self.OutPath + fileName + ".dat")            

        nprint(" done.\n")
        nprint(" ------------------------------------------------------------------")
    
    def run(self):
        """Runs the sumulation."""

        self.clock.tic()
        
        nprint("\n\n Starting simulation ...")

        # Print simulation info
        nprint(self.staticnstr())
        eprint(self.debugstr())
             
        if (self.check() == False):
            raise error(" ERROR: Check failed !!!")        
        
        # Debug print
        dprint("\n" + str(self.VDD))
        dprint("\n" + str(self.VGG))

        # save simulation parameters
        if (self.workers.IAmMaster()):
            # Create directory if not exist.
            if not os.path.exists(self.OutPath):
                os.makedirs(self.OutPath)            
            # save simulation parameters to a pickle file
            fp = open(self.OutPath + self.OutFileName + ".pkl", 'wb')
            pk.dump(self, fp)
            fp.close()                        

        # Save electrostatic geometry
        if (self.verbosity >= vprint.MSG_DEBUG):
            if (self.workers.IAmMaster()):
                # Create directory if not exist.
                if not os.path.exists(self.OutPath):
                    os.makedirs(self.OutPath) 
                self.V.exportSvg(self.OutPath + self.DebugGeomFile)
        
        # Configure the RGF solver
        if self.kp.N == 0: # without k-loop
            self.rgf = CohRgfLoop(self.workers, self.nb, self.kT, self.ieta,  
                    self.OrthoBasis)
        else:# with k-loop
            self.rgf = CohRgfLoop(self.workers, self.nb, self.kT, self.ieta,  
                    self.OrthoBasis, 2)

        # Setup H and S 
        if (self.DevType == self.COH_RGF_UNI):        # for uniform RGF blocks

            for ib in range(0, self.nb+1):            # setup the block hamiltonian
                if self.kp.N == 0: # no k-loop
                    if (ib != self.nb):
                        self.rgf.H0(self.H0, ib)          # H0: 0 to N+1
                        self.rgf.S0(self.S0, ib)          # S0: 0 to N+1
                    self.rgf.Hl(self.Hl, ib)              # Hl: 0 to N+2
                    
                else: # we have k-loop, add transverse neighbors
                    self.rgf.k(self.kp.kp)                      # set k-points
                    if (ib != self.nb):
                        self.rgf.H0(self.H0[0], ib, 0)          # H0_i,i: 0 to N+1 
                        self.rgf.H0(self.H0[1], ib, 1)          # H0_i,i+1: 0 to N+1
                        self.rgf.H0(self.H0[2], ib, 2)          # H0_i,i-1: 0 to N+1
                        self.rgf.S0(self.S0[0], ib, 0)          # S0_i,i: 0 to N+1
                        self.rgf.pv0(self.pv[0], ib, 0)
                        self.rgf.pv0(self.pv[1], ib, 1)
                        self.rgf.pv0(self.pv[2], ib, 2)
                    self.rgf.Hl(self.Hl[0], ib, 0)              # Hl_i,i: 0 to N+2
                    self.rgf.Hl(self.Hl[1], ib, 1)              # Hl_i,i: 0 to N+2
                    self.rgf.Hl(self.Hl[2], ib, 2)              # Hl_i,i: 0 to N+2
                    self.rgf.pvl(self.pvl[0], ib, 0)
                    self.rgf.pvl(self.pvl[1], ib, 1)
                    self.rgf.pvl(self.pvl[2], ib, 2)                
                    
        elif (self.DevType == self.COH_RGF_NON_UNI):  # Non-uniform RGF blocks        
            for ib in range(0, self.nb):                 # setup the block hamiltonian
                self.rgf.H0(self.H0[ib], ib)             # H0: 0 to N+1=nb-1
                self.rgf.S0(self.S0, ib)              # S0: just the identity matrix stored in S0(0)
                if (ib > 0):
                    self.rgf.Hl(self.Hl[ib], ib)     # Hl: 1 to N+1=nb-1                            
            self.rgf.Hl(self.Hl[1], 0)               # Set H_0,-1 = H_1,0. Hl(0) = H_0,-1
            self.rgf.Hl(self.Hl[ib], ib+1)           # Set H_N+2,N+1 = H_N+1,N            
        # Enable calculations
        for type, value in self.Calculations.iteritems():
            # transmission
            if (type == "TE"):
                self.rgf.enableTE(value)
            # current
            if (type == "I"):
                if (isinstance( value, int)):
                    self.rgf.enableI(value, 0, 1)
                else:
                    for I in self.Calculations["I"]:
                        self.rgf.enableI(I["N"], I["From"], I["To"])
            if (type == "DOS"):
                self.rgf.enableDOS(value)
            if (type == "n"):
                if (isinstance( value, int)):
                    self.rgf.enablen(value)
                else:
                    for n in self.Calculations["n"]:
                        if (n["Block"] == "All"):
                            for ib in range(1, self.nb-2):
                                self.rgf.enablen(n["N"], ib)
                        else:
                            self.rgf.enablen(n["N"], n["Block"])
    
        # Loop over drain and gate bias
        for VDD in self.VDD:
            for VGG in self.VGG:
                self.runBiasStep(VGG, self.Vo, VDD)
                pass
        self.clock.toc()           
        nprint("\n" + str(self.clock) + "\n")
 
    def updateBoundingBox(self):
        """Updates the bounding box of atomistic geometry."""
        # Just to make sure that no point of gate regions is    
        # at the border
        a = self.hp.a
        delta = a*7.0/220.0
        self.xmn = self.geom.xmin - delta   
        self.xmx = self.geom.xmax + delta
        self.ymn = self.geom.ymin - delta
        self.ymx = self.geom.ymax + delta        
        
    def __str__(self):
        msg = self.staticnstr()
        msg += self.dynamicnstr()
        return msg
    
    def staticnstr(self):
        msg = "\n Transport simulation parameters:"
        
        msg += "\n  Device geometry:"
        msg += "\n  ---------------------------------------------"
        msg += "\n  ... |-1 | 0 | 1 |   ...  | " + str(self.nb-2) + " | " + str(self.nb-1) + " | " + str(self.nb) + " | ..."
        msg += "\n  ---------------------------------------------"
        msg += "\n            ^  <-------^------->  ^"
        msg += "\n          left       Device     right"
        msg += "\n        contact                contact"
        
        # Bias information
        msg += "\n Bias: "
        msg += "\n  VDD: min " + str(min(self.VDD)) + ", max " + str(max(self.VDD)) 
        msg += ", number " + str(len(self.VDD))
        msg += "\n  VGG: min " + str(min(self.VGG)) + ", max " + str(max(self.VGG)) 
        msg += ", number " + str(len(self.VGG))
        msg += "\n  Vo: " + str(self.Vo)
        
        # Device information
        msg += "\n Device:"
        msg += "\n  Lengh: " + str(self.geom.xl) + "(" + str(self.geom.xl/10.0) + " nm)"
        msg += ", Width: " +  str(self.geom.yl) + "(" + str(self.geom.yl/10.0) + " nm)"
        msg += ", Height: " + str(self.geom.zl)+ "(" + str(self.geom.zl/10.0) + " nm)"
        msg += "\n  Num of atoms: " + str(self.geom.NumOfAtoms) 
        msg += ", Num of orbitals: " +  str(self.geom.NumOfOrbitals)
        
        # Print Hamiltonian parameters.
        msg += "\n"
        msg += str(self.hp)
       
        # Calculations to perform
        msg += " Calculations:\n"
        for type, value in self.Calculations.iteritems():
            if (type == "TE"):
                msg += "  Transmission.\n"
            if (type == "I"):
                for I in self.Calculations["I"]:
                    msg += "  Current from block # " + str(I["From"])
                    msg += " to block # " + str(I["To"]) 
                    msg += " (" + str(I["N"]) + ").\n"
            if (type == "n"):
                msg += "  Electron density.\n"
                
        msg += "  Save output at: " + self.OutPath + self.OutFileName + "*\n"
                    
        # End of info section
        msg += " ------------------------------------------------------------------"
        
        return msg
    
    def dynamicnstr(self):
        # NEGF parameters
        msg = str(self.rgf)
        # Electrostatics
        msg += "\n"
        msg += str(self.V)

        return msg
 

    def debugstr(self):
        msg = " Debugging information: \n"
        msg = "  Atomic structure: \n"
        msg += str(self.geom)
        msg += " ------------------------------------------------------------------"
        return msg
        
    def __getstate__(self):
        dct = dict(self.__dict__)
        del dct['workers']
        del dct['clock']
        del dct['kp']
        #del dct['rgf']
        return dct
    
    def __setstate__(self, dct):
        self.__dict__.update(dct)
        
"""
 Loads and returns a Transport object from pickle file
"""
def loadTransport(fileName):
    pf = open(fileName)
    tr = pk.load(pf)
    pf.close()
    return tr






