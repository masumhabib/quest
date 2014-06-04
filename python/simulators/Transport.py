"""
    Qunantum Transport simulator using QMICAD.
    Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>

    Last update: 05/12/2014
"""

import os
import pickle as pk
import numpy as np
import random as rn

#import qmicad
from qmicad.atoms import AtomicStruct, SVec
from qmicad.hamiltonian import TISurfKpParams4, TISurfKpHam4, TISurfKpParams, TISurfKpHam, GrapheneKpParams, GrapheneKpHam
from qmicad.negf import CohRgfaParams, NegfEloop
from qmicad.potential import LinearPot
from qmicad.utils import VecGrid, Timer, Workers, vprint
from qmicad._utils.vprint import nprint, dprint, eprint

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

        self.version = qmicad.version # Library version
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
        self.nbw        = []           # Width of all layers

        # Hamiltonian type
        self.HAM_TI_SURF_KP  = 10       # TI surface k.p hamiltonian
        self.HAM_TI_SURF_KP4 = 11       # TI surface k.p hamiltonian with 4 spin basis set
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
        
        # Debug stuffs
        self.DebugPotFile   = "dbg_pot.dat"
        self.DebugGeomFile  = "dbg_geom.svg"
        
        # Dry run
        self.DryRun         = False
        
        # Skip if resulting dat file already exists?
        self.SkipExistingSimulation = False

    @property 
    def verbosity(self):
        return self._verbosity  
    @verbosity.setter
    def verbosity(self, value):
        self._verbosity = value
        # Set verbosity level of QMICAD library.
        qmicad.setVerbosity(self._verbosity)
        vprint.verbosity = self._verbosity
   
    def muS(self, VDD):
        """Calculates the fermi energy level of the source."""
        return self.np.mu - VDD*self.V.rVS
    
    def muD(self, VDD):
        """Calculates the drain energy level of the source."""
        return self.np.mu - VDD*self.V.rVD
 
    
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
        
        # TI k.p surface      
        if (self.HamType == self.HAM_TI_SURF_KP):    
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = TISurfKpParams()
            else:
                if not isinstance(self.hp, TISurfKpParams): # if hp exists but not TISurfKpParams type, create it.
                    self.hp = TISurfKpParams()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct()
            self.geom.genRectLattAtoms(self.nb, self.nw, self.hp.a, self.hp.a, self.hp.ptable)
            # Just to make sure that no point of gate regions is    
            # at the border
            a = self.hp.a
            delta = a*7.0/220.0
            self.nbw = [self.nw]*self.nb
        elif (self.HamType == self.HAM_TI_SURF_KP4):    
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = TISurfKpParams4()
            else:
                if not isinstance(self.hp, TISurfKpParams4): # if hp exists but not TISurfKpParams type, create it.
                    self.hp = TISurfKpParams4()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct()
            self.geom.genRectLattAtoms(self.nb, self.nw, self.hp.a, self.hp.a, self.hp.ptable)
            # Just to make sure that no point of gate regions is    
            # at the border
            a = self.hp.a
            delta = a*7.0/220.0  + a*7.0/2200.0
            self.nbw = [self.nw]*self.nb
        elif (self.HamType == self.HAM_GRAPHENE_KP):
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = GrapheneKpParams()
            else:
                if not isinstance(self.hp, GrapheneKpParams): # if hp exists but not GrapheneKpParams type, create it.
                    self.hp = GrapheneKpParams()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct()
            self.geom.genRectLattAtoms(self.nb, self.nw, self.hp.a, self.hp.a, self.hp.ptable)
            # Just to make sure that no point of gate regions is    
            # at the border
            a = self.hp.a
            delta = a*7.0/220.0
            self.nbw = [self.nw]*self.nb
        else:
            raise RuntimeError(" Unsupported Hamiltonian type. ")
        
        self.updateBoundingBox()

    def createRoughEdges(self, sigma):
        """ 
        Creates rough edges. Works only for sorted lattice points.
        For rectangular lattice 
        """
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
                
    def generateHamiltonian(self):
        """ Generates hamiltonian and overlap matrices. """
        # Create Hamiltonian generator
        if (self.HamType == self.HAM_TI_SURF_KP):
            self.ham = TISurfKpHam(self.hp)   
        if (self.HamType == self.HAM_TI_SURF_KP4):
            self.ham = TISurfKpHam4(self.hp)               
        if (self.HamType == self.HAM_GRAPHENE_KP):
            self.ham = GrapheneKpHam(self.hp)   

        # For uniform RGF blocks
        if (self.DevType == self.COH_RGF_UNI):            
            lyr0 = self.geom.span(0, self.nw-1);               # extract block # 0
            lyr1 = self.geom.span(self.nw, 2*self.nw-1);    # extract block # 1
            # Block hamiltonian.
            self.ham.setSizeForNegf(2)                            # uniform device, only need H_0,0 and H_1,0
            self.ham.genDiagBlock(lyr0, lyr0, 0)                  # generate H_0,0 and S_0,0
            self.ham.genLowDiagBlock(lyr1, lyr0, 0)               # generate H_1,0                 
            # Set block Hamiltonian for the NEGF calculation.
            self.np = CohRgfaParams(self.nb)            
            for ib in range(0, self.nb+1):                  # setup the block hamiltonian
                if (ib != self.nb):
                    self.np.H0(self.ham.H0(0), ib)          # H0: 0 to N+1
                    self.np.S0(self.ham.S0(0), ib)          # S0: 0 to N+1
                self.np.Hl(self.ham.Hl(0), ib)              # Hl: 0 to N+2
        # Non-uniform RGF blocks        
        elif (self.DevType == self.COH_RGF_NON_UNI):
            # Set block Hamiltonian for the NEGF calculation.
            self.ham.setSizeForNegf(self.nb)            # non-uniform device, store all H_i,i and H_i,i-1
            self.np = CohRgfaParams(self.nb)            
            beg = 0
            for ib in range(0, self.nb):                    # setup the block hamiltonian
                end = beg + self.nbw[ib] - 1
                
                # generate H_i,i and S_i,i
                lyri = self.geom.span(beg, end);            # extract block # i             
                self.ham.genDiagBlock(lyri, lyri, ib)
                self.np.H0(self.ham.H0(ib), ib)             # H0: 0 to N+1=nb-1
                self.np.S0(self.ham.S0(0), ib)              # S0: just the identity matrix stored in S0(0)

                # generate H_i,i-1
                if (ib > 0):
                    self.ham.genLowDiagBlock(lyri, lyrim1, ib)      
                    self.np.Hl(self.ham.Hl(ib), ib)              # Hl: 1 to N+1=nb-1                
                
                lyrim1 = lyri
                beg = end + 1
            
            self.np.Hl(self.ham.Hl(1), 0)               # Set H_0,-1 = H_1,0. Hl(0) = H_0,-1
            self.np.Hl(self.ham.Hl(ib), ib+1)           # Set H_N+2,N+1 = H_N+1,N

        
    def setupPotential(self):
        """ Sets up the potential profile """
        # Linear
        if (self.PotType == self.POT_LINEAR):
            self.V = LinearPot(self.geom)
    
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
        self.np.muS = self.muS(VDD)
        self.np.muD = self.muD(VDD)
            
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
#        if (self.DevType == self.COH_RGF_UNI): 
        beg = 0
        for ib in range(self.nb):                  # setup the block hamiltonian
            end = beg + self.nbw[ib] - 1
            self.np.V(self.V.toOrbPot(beg, end), ib)
            beg = end + 1
            
#                self.np.V(self.V.toOrbPot(self.nw*ib, 
#                                          self.nw*(ib+1)-1), ib)

        # Create energy grid
        if (self.np.AutoGenE):
            Emin = self.np.muD - 10*self.np.kT
            Emax = self.np.muS + 10*self.np.kT
        else:
            Emin = self.np.Emin
            Emax = self.np.Emax
            
        EE = VecGrid(Emin, Emax, self.np.dE)

        
        nprint("\n Energy grid:\n  Min: " + str(Emin) + ", max: " 
                    + str(Emax) + ", interval " + str(self.np.dE)
                    + ".")
        nprint("\n Total " + str(EE.N()) + " energy point(s) " 
                    + "running on " + str(self.workers.N()) + " CPU(s): " 
                    + str(int(round(EE.N()/self.workers.N()))) + " pts/CPU ... \n")
                                
        # create and run energy loop
        Eloop = NegfEloop(EE, self.np, self.workers) 

        # enable calculations
        for type, value in self.Calculations.iteritems():
            # transmission
            if (type == "TE"):
                Eloop.enableTE(value)
            # current
            if (type == "I"):
                for I in self.Calculations["I"]:
                    Eloop.enableI(I["N"], I["From"], I["To"])
            if (type == "DOS"):
                Eloop.enableDOS(value)
            if (type == "n"):
                Eloop.enablen(value)
        
        # Run the simulation
        if (self.DryRun == False):
            Eloop.run()
        nprint("\n  done.")
        
        nprint("\n Saving results to disk ...")
        # save results
        if (self.workers.IAmMaster()):
            # Create directory if not exist.
            if not os.path.exists(self.OutPath):
                os.makedirs(self.OutPath)            
            # save results to file.
            fo = open(self.OutPath + fileName + ".dat", "wt")
            fo.close()
            Eloop.save(self.OutPath + fileName + ".dat")            

        nprint(" done.\n")
        nprint(" ------------------------------------------------------------------")
    
    def run(self):
        """Runs the sumulation."""

        self.clock.tic()
        
        # Print welcome message.
        nprint(qmicad.greet())
        
        # Print simulation info
        nprint(self.nstr())
        eprint(self.dstr())
             
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
        
        # Loop over drain bias
        for VDD in self.VDD:
            for VGG in self.VGG:
                self.runBiasStep(VGG, self.Vo, VDD)
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
        return self.nstr()
    
    def nstr(self):
        msg = "\n Transport simulation parameters:"
        
        msg += "\n  Device geometry:"
        msg += "\n  ---------------------------------------------"
        msg += "\n  ... |-1 | 0 | 1 |   ...  | " + str(self.np.nb-2) + " | " + str(self.np.nb-1) + " | " + str(self.np.nb) + " | ..."
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
        
        # NEGF parameters
        msg += str(self.np)
        
        # Electrostatics
        msg += "\n"
        msg += str(self.V)
        
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
        msg += "  Save output at: " + self.OutPath + self.OutFileName + "*\n"
                    
        # End of info section
        msg += " ------------------------------------------------------------------"
        
        return msg
    
    def dstr(self):
        msg = " Debugging information: \n"
        msg = "  Atomic structure: \n"
        msg += str(self.geom)
        msg += " ------------------------------------------------------------------"
        return msg
        
    def __getstate__(self):
        dct = dict(self.__dict__)
        del dct['workers']
        del dct['clock']
        del dct['ham']
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






