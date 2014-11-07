
"""
    Band structure calculator using QMICAD.
    Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>

    Last update: 05/07/2014
"""

import os
import pickle as pk
import numpy as np
from   math import pi

from qmicad import setVerbosity, greet
from qmicad.atoms import AtomicStruct, LCoord
from qmicad.hamiltonian import TISurfKpParams4, TISurfKpParams, TI3DKpParams, GrapheneKpParams, GrapheneTbParams, generateHamOvl
from qmicad.kpoints import KPoints
from qmicad.band import BandStruct
from qmicad.utils import Timer, Workers, Point
from qmicad.simulators.utils import vprint
from qmicad.simulators.utils.vprint import nprint, dprint, eprint


class Band(object):
    """
    Band structure calculator using QMICAD.    
    """
    
    def __init__(self, workers = None):   
        """ 
        Class constructor. 
        param:
          workers: mpi.world - MPI world communicator.
        """  

        self.verbosity = vprint.MSG_NORMAL # Verbosity level
        
        # MPI stuff
        if workers is None:
            self.workers = Workers()
        else:
            self.workers = Workers(workers)
        vprint.IAmMaster = self.workers.IAmMaster()
    
        # Timer
        self.clock = Timer()

        # output path settings
        self.OutPath    = "./out/"     # Output path
        self.OutFileName = "EK"        # Output file name prefix
                        
        # Cell parameters
        self.nl         = 1            # Length of the cell
        self.nw         = 1            # Width of the cell
        self.nh         = 1            # Height of the cell
        
        # Band structure simulation parameters
        self.Dim        = 1            # Cell dimension
        self.nb         = 2            # Number of bands to be saved
        self.nk1        = 100          # number of k points along a1
        self.nk2        = 100          # number of k points along a2
        
        self.lc         = []           # Neighbor list
        self.kp         = KPoints()    # K point generator

        # Hamiltonian type
        self.HAM_TI_SURF_KP  = 10       # TI surface k.p hamiltonian
        self.HAM_TI_SURF_KP4 = 11       # TI surface k.p hamiltonian with 4 spin basis set
        self.HAM_TI_3D_KP    = 15       # TI 3D k.p hamiltonian with 4 spin basis set
        self.HAM_GRAPHENE_KP = 20       # Graphene k.p Hamiltonian
        self.HamType         = self.HAM_TI_SURF_KP 
                        
        # Calculations
        self.EnableEigVec = False
                
        # Dry run
        self.DryRun             = False
        
        # Skip if resulting dat file already exists?
        self.SkipExistingSimulation = False        

    @property 
    def verbosity(self):
        return self._verbosity  
    @verbosity.setter
    def verbosity(self, value):
        self._verbosity = value
        # Set verbosity level of QMICAD library.
        setVerbosity(self._verbosity)
        vprint.verbosity = self._verbosity

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
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nl, self.nw, self.nh)
        elif (self.HamType == self.HAM_TI_SURF_KP4):    
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = TISurfKpParams4()
            else:
                if not isinstance(self.hp, TISurfKpParams4): # if hp exists but not TISurfKpParams type, create it.
                    self.hp = TISurfKpParams4()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nl, self.nw, self.nh)
        elif (self.HamType == self.HAM_TI_3D_KP):    
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = TI3DKpParams()
            else:
                if not isinstance(self.hp, TI3DKpParams): # if hp exists but not TISurfKpParams type, create it.
                    self.hp = TI3DKpParams()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nl, self.nw, self.nh)
        elif (self.HamType == self.HAM_GRAPHENE_KP):
            # Hamiltonian parameter
            if not hasattr(self, "hp"): # if hp does not exist, create it.
                self.hp = GrapheneKpParams()
            else:
                if not isinstance(self.hp, GrapheneKpParams): # if hp exists but not GrapheneKpParams type, create it.
                    self.hp = GrapheneKpParams()
            # Create atomistic geometry of the device.
            self.geom = AtomicStruct(self.hp.ptable)
            self.geom.genSimpleCubicStruct(self.hp.ptable[0], self.hp.a, self.nl, self.nw, self.nh)
        else:
            raise RuntimeError(" Unsupported Hamiltonian type. ")
        
        # Save lattice vector
        self.lv = self.geom.LatticeVector                  

    def addNeighbor(self, n1 = 0, n2 = 0, n3 = 0, AutoGen = True):
        if AutoGen == True:
            if (self.Dim == 1):
                # The nearest neighbor on the right,
                # need H_0,0 and H_0,1 only.
                self.lc.append(LCoord(0, 0, 0))
                self.lc.append(LCoord(1, 0, 0))

            elif (self.Dim == 2):
                # For rectangular lattice.
                if (self.HamType == self.HAM_TI_SURF_KP or 
                    self.HamType == self.HAM_TI_SURF_KP4 or 
                    self.HamType == self.HAM_TI_3D_KP or 
                    self.HamType == self.HAM_GRAPHENE_KP):
                    # The nearest neighbors
                    self.lc.append(LCoord(0, 0, 0))                     # main cell
                    self.lc.append(LCoord(1, 0, 0))                     # neighbor on the right      
                    self.lc.append(LCoord(1, 1, 0))                     # neighbor on the right-top      
                    self.lc.append(LCoord(0, 1, 0))                     # neighbor on the top      
                    self.lc.append(LCoord(-1, 1, 0))                    # neighbor on the top-left                       
                else:
                    raise RuntimeError(" Unsupported lattice type. ")
            
        else:
            # Add a neighbor in terms of lattice vector
            self.lc.append(lcoord(n1, n2, n3))
        
    def generateHamiltonian(self):
        """ Generates hamiltonian and overlap matrices. """

        # Generate hamiltonian and overlap matrices.
        # Hamiltonian of nearest neighbors
        nn = len(self.lc)                                 # number of nearest neighbors.
        # generate H_0,i and S_0,i
        self.H = [];
        self.S = [];
        all_neigh = AtomicStruct()                       # for storing all neighbors
        for inn in range(nn):
            neigh = self.geom + self.lc[inn]                  # extract block # 1
            H, S = generateHamOvl(self.hp, self.geom, neigh)
            self.H.append(H)
            all_neigh = all_neigh + neigh                 # collect neighbors for debug.            
        if self.verbosity == vprint.MSG_DUMP:
            if (self.workers.IAmMaster()):
                if not os.path.exists(self.OutPath):
                    os.makedirs(self.OutPath)    
                nprint (" Saving atomistic geometry to " + self.OutPath + "dbg_geom.gjf")
                all_neigh.exportGjf(self.OutPath + 'dbg_geom.gjf')
      
    def generateKPoints(self):
        if (self.Dim == 1):
            la1 = self.lv.la1()
            self.kp.addKLine(Point(-pi/la1, 0), Point(pi/la1, 0), self.nk1)
        elif (self.Dim == 2):
            # For rectangular lattice.
            if (self.HamType == self.HAM_TI_SURF_KP or 
                self.HamType == self.HAM_TI_SURF_KP4 or 
                self.HamType == self.HAM_GRAPHENE_KP):
                la1 = self.lv.la1()
                la2 = self.lv.la2();
                self.kp.addKRect(Point(-pi/la1, -pi/la2), Point(pi/la1, pi/la2), self.nk1, self.nk2)                
            else:
                raise RuntimeError(" Unsupported lattice type. ")

    def check(self):
        ret = True            
        return ret
        
    def run(self):
        """ Runs the sumulation. """
        self.clock.tic()
        
        # Print welcome message.
        nprint(greet())
        
        # Print simulation info
        nprint(self.nstr())
        eprint(self.dstr())
             
        if (self.check() == False):
            raise error(" ERROR: Check failed !!!")        
        
        # save simulation parameters
        if (self.workers.IAmMaster()):
            # Create directory if not exist.
            if not os.path.exists(self.OutPath):
                os.makedirs(self.OutPath)            
            # save simulation parameters to a pickle file
            fp = open(self.OutPath + self.OutFileName + ".pkl", 'wb')
#            pk.dump(self, fp)
            fp.close()
            
        # Setup band structure calculator.
        nn = len(self.lc)
        bs = BandStruct(self.workers, nn)
        bs.lv = self.lv                            # Lattice vector.
        bs.nb = self.nb                            # number of bands
        bs.ne = self.geom.NumOfElectrons           # number of bands        
        for inn in range(nn):
            bs.lc(self.lc[inn], inn)               # lattice coordinate
            bs.H(self.H[inn], inn)                 # Hamiltonian
        bs.k(self.kp.kp)                           # kpoints

        # Calculate eigen vectors?
        if self.EnableEigVec == True:
            bs.enableEigVec()
                
        nprint("\n Total " + str(self.kp.N) + " k point(s) " 
            + "running on " + str(self.workers.N()) + " CPU(s)...\n")

        # skip calculation if result file exists.
        if self.SkipExistingSimulation == True:
            fileName = self.OutFileName + ".dat"
            if os.path.isfile(self.OutPath + fileName):
                nprint("\n  Result exists, skipping.")
                return
            
        # Run the simulation
        if (self.DryRun == False):
            bs.run()
        nprint("\n  done.")
         
        nprint("\n Saving results to disk ...")
        # save results
        if (self.workers.IAmMaster()):
            # Create directory if not exist.
            if not os.path.exists(self.OutPath):
                os.makedirs(self.OutPath)
            # save to file.
            fileName = self.OutFileName + ".dat"
            bs.save(self.OutPath + fileName)
             
        nprint(" done.")
        nprint("\n ------------------------------------------------------------------")
        
        self.clock.toc()           
        nprint("\n" + str(self.clock) + "\n")

    def __str__(self):
        return self.nstr()
    
    def nstr(self):
        msg = "\n Band structure simulation parameters:"
                
        # Unit cell information
        msg += "\n Cell:"
        msg += "\n  Lengh: " + str(self.geom.xl) + "(" + str(self.geom.xl/10.0) + " nm)"
        msg += ", Width: " +  str(self.geom.yl) + "(" + str(self.geom.yl/10.0) + " nm)"
        msg += ", Height: " + str(self.geom.zl)+ "(" + str(self.geom.zl/10.0) + " nm)"
        msg += "\n  Num of atoms: " + str(self.geom.NumOfAtoms) 
        msg += ", Num of orbitals: " +  str(self.geom.NumOfOrbitals)
        msg += "\n  Num of electrons: " + str(self.geom.NumOfElectrons) 
        
        # Print Hamiltonian parameters.
        msg += "\n"
        msg += str(self.hp)
        
        # Band structure parameters
#        msg += str(self.bp)
                
        # Calculations to perform
        msg += " Calculations:\n"
        msg += "  EK.\n"
        if self.EnableEigVec == True:
            msg += "  Eigen Vectors.\n"

        msg += "  Save output at: " + self.OutPath + self.OutFileName + "*\n"
                    
        # End of info section
        msg += " ------------------------------------------------------------------"
        
        return msg
    
    def dstr(self):
        msg = "\n Debugging information: \n"
        msg += str(self.geom)
        msg += " ------------------------------------------------------------------"
        return msg
        
    def __getstate__(self):
        dct = dict(self.__dict__)
        del dct['workers']
        del dct['clock']
        del dct['ham']
        del dct['kp']
        del dct['bp']
        return dct
    
    def __setstate__(self, dct):
        self.__dict__.update(dct)
        
"""
 Loads and returns a Band object from pickle file
"""
def loadBand(fileName):
    pf = open(fileName)
    bd = pk.load(pf)
    pf.close()
    return bd










