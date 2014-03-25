"""
    Qunantum Transport simulator using QMICAD.
     
    Author: K M Masum Habib <masum.habib@virginia.edu>
    Last update: 02/14/2014
"""

import os
import pickle as pk
import numpy as np

import qmicad
from qmicad.atoms import AtomicStruct
from qmicad.hamiltonian import TISurfKpParams, TISurfKpHam
from qmicad.negf import CohRgfaParams, NegfEloop
from qmicad.potential import LinearPot
from qmicad.utils import VecGrid, Timer, Workers, vprint
from qmicad.utils.vprint import nprint, dprint, eprint

"""
    Generic transport class.
"""
class Transport(object):
    ## 
    # class initialization. 
    # @Parameters:
    #    sp: TiSimuParams - TI simulation parameters
    #    
    def __init__(self, workers):   
        # Library version
        self.version = qmicad.version
        
        # Verbosity level
        self.verbosity = vprint.MSG_NORMAL

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
        self.nw         = 5            # Width of the device

        # Hamiltonian type
        self.HAM_TI_SURF_KP = 10       # TI surface k.p hamiltonian
        self.HamType        = self.HAM_TI_SURF_KP 
        
        # Device types
        self.COH_RGF_UNI    = 10       # Coherent Uniform RGF type device
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

    @property 
    def verbosity(self):
        return self._verbosity  
    @verbosity.setter
    def verbosity(self, value):
        self._verbosity = value
        # Set verbosity level of QMICAD library.
        qmicad.setVerbosity(self._verbosity)
        vprint.verbosity = self._verbosity
   
    # Source fermi energy
    def muS(self, VDD):
        return self.np.mu - VDD*self.V.rVS
    
    # Drain fermi energy
    def muD(self, VDD):
        return self.np.mu - VDD*self.V.rVD
 
    # Gate voltage of gate #ig
    def VG(self, VGG, Vo, ig):
        # Set gate voltages
        return VGG + Vo*self.V.rVG[ig]
      
    def addSource(self, sql, rVS = -0.5):
        self.V.addSource(sql)
        self.V.rVS = rVS
        # just for pickling
        self.V.sql = sql    
    
    def addDrain(self, dql, rVD = 0.5):
        self.V.addDrain(dql)
        self.V.rVD = rVD
        # just for pickling
        self.V.dql = qql 

    def addGate(self, gql, rVG = 1):
        self.V.addGate(gql)
        self.V.rVG.append(rVG)
         # just for pickling
        self.V.gql.append(gql) 

    def addLinearRegion(self, lql):
        self.V.addLinearRegion(lql)
        # just for pickling
        self.V.lql.append(lql) 
    
    # Atomistic geometry            
    def createAtomicGeom(self):
        # TI k.p surface      
        if (self.HamType == self.HAM_TI_SURF_KP):    
            # Hamiltonian parameter
            self.hp = TISurfKpParams()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct(self.nb, self.nw, self.hp.ax, self.hp.ay, self.hp.ptable)

    # Generate hamiltonian and overlap matrices.
    def generateHamiltonian(self):
        # Create Hamiltonian generator
        self.ham = TISurfKpHam(self.hp)        

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

    # Potential profile        
    def setupPotential(self):
        # Linear
        if (self.PotType == self.POT_LINEAR):
            self.V = LinearPot(self.geom)


    # Prepare simulation
    def prepare(self):
        self.createAtomicGeom()
        self.generateHamiltonian()
        self.setupPotential()

    # 
    def check(self):
        return True 

    ## 
    # Runs the sumulation.
    #
    def runBiasStep(self, VGG, Vo, VDD):            
        # Set drain and Fermi levels
        self.np.muS = self.muS(VDD)
        self.np.muD = self.muD(VDD)
            
        nprint("\n Bias loop:")                                
        nprint("\n  VGG = " + str(VGG) + ", Vo = " + str(Vo) 
                    + ", VDD = " + str(VDD) + ".")
        
        fileName = self.OutFileName + "_VGG{0:2.3f}_Vo{1:2.3f}_VDD{2:2.3f}".format(VGG, Vo, VDD)

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
                if (ig == len(self.rVG) - 1):
                    self.V.VD(VG)
        
        # Potential of linear region
        if (self.PotType == self.POT_LINEAR):
            for il in range(self.V.NLR):
                VG1 = self.VG(VGG, il)
                VG2 = self.VG(VGG, il+1)
                self.V.VLR(il, VG1, VG2)        

        # Compute the electrostatic potential
        self.V.compute()
        eprint(self.V)
        if (self.verbosity >= vprint.MSG_DEBUG):
            self.V.exportPotential(self.OutPath + self.DebugPotFile)
        
        # Export potential to NEGF
        if (self.DevType == self.COH_RGF_UNI): 
            for ib in range(self.nb):                  # setup the block hamiltonian
                self.np.V(self.V.toOrbPot(self.nw*ib, 
                                          self.nw*(ib+1)-1), ib)

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
                    + str(EE.N()/self.workers.N()) + " pts/CPU ... \n")
                        
        # create and run energy loop
        Eloop = NegfEloop(EE, self.np, self.workers) 

        # enable calculations
        for type, value in self.Calculations.iteritems():
            # transmission
            if (type == "TE"):
                Eloop.enableTE(value)
            # current
            if (type == "I"):
                for ib, N in self.Calculations["I"].iteritems():
                    Eloop.enableI(ib, N)
            if (type == "DOS"):
                Eloop.enableDOS(value)
            if (type == "n"):
                Eloop.enablen(value)
        
        # Run the simulation
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
            
            # save simulation parameters to a pickle file
            fp = open(self.OutPath + fileName + ".pkl", 'wt')
            pk.dump(self, fp)
            fp.close()

        nprint(" done.\n")
        nprint(" ------------------------------------------------------------------")
    
    ## 
    # Runs the sumulation.
    #
    def run(self):
        self.clock.tic()
        
        # Print welcome message.
        nprint(qmicad.greet())
        nprint(self.__str__())
             
        if (self.check() == False):
            raise error(" ERROR: Check failed !!!")        

        if (self.verbosity >= vprint.MSG_DEBUG):
            self.V.exportSvg(self.sp.OutPath + self.DebugGeomFile)
        
        # Debug print
        dprint(self.VDD)
        dprint(self.VGG)

        # Loop over drain bias
        for VDD in self.VDD:
            for VGG in self.VGG:
                self.runBiasStep(VGG, self.Vo, VDD)
        self.clock.toc()           
        nprint("\n" + str(self.clock) + "\n")
 
    def __str__(self):
        msg = "\n Transport simulation parameters:"
        
        # Bias information
        msg += "\n Bias: "
        msg += "\n  VDD: min " + str(min(self.VDD)) + ", max " + str(max(self.VDD)) 
        msg += ", number " + str(len(self.VDD))
        msg += "\n  VGG: min " + str(min(self.VGG)) + ", max " + str(max(self.VGG)) 
        msg += ", number " + str(len(self.VGG))
        msg += "\n  Vo: " + str(self.Vo)
        
        # Device information
        msg += "\n Device:"
        msg += "\n  Lengh: " + str(self.geom.xl) + ", Width: " +  str(self.geom.yl)
        msg += ", Height: " + str(self.geom.zl)
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
                for ib, N in self.Calculations["I"].iteritems():
                    msg += "  Current at block #" + str(ib) + ".\n"
        msg += "  Save output at: " + self.OutPath + self.OutFileName + "*\n"
                    
        # End of info section
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
def load(fileName):
    pf = open(fileName)
    tr = pickle.load(pf)
    pf.close()
    return tr






