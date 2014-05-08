## 
# Topological Insulator surface band structure simulation parameters: 
# contains all the simulation parameters needed for k.p
# and band structure model.
# Instance of thid class is required by all other classes.
#
# Author: K M Masum Habib <masum.habib@virginia.edu>
# Last Update: 02/14/2014
#

import os
import  pickle
from    qmicad     import *
from    math       import *

class TiSurfBandSimuParams:
    
    ## 
    # class initialization. 
    # @Parameters:
    #    workers: mpi.world - MPI world communication between worker processes.
    #    
    def __init__(self, workers):         
        
        
        # output settings
        self.OutPath     = "./band/"    # Output path
        self.OutFileName = "ti"        # Output file prefix
        
        # Device parameters
        self.nl         = 1            # Length of the device+contacts
        self.nw         = 3            # Width of the device
        
        # Band structure simulation parameters
        self.Dim        = 1            # Dimension: 1D/2D?
        self.mu         = 0.0          # Device Fermi level
        self.nb         = 2            # Number of bands to be saved
        self.dk         = 0.01         # del k

        # TI surface k.p parameters
        self.a          = 2.0          # k.p grid spacing
        self.K          = 0.57971      # Coefficient ot solve Fermion doubling
        self.dtol       = 1E-3         # Tolerance when considering neighbors
        self.A2         = 3.33         # A2 paramter
        self.C          = 0.0          # C parameter

        # Periodic table with fake atoms
        self.ptable     = PeriodicTable()       # Fake periodic table for k.p
        self.ptable.add(0, "D", 2, 2)           # Fake atom
        
        # Verbosity level
        self.MSG_QUIET  = 0             # Do not print any message
        self.MSG_NORMAL = 1             # Print some messages
        self.MSG_DEBUG  = 10            # Print debug message
        self.MSG_DUMP   = 100           # Print everyting
        self.verbosity    = self.MSG_DEBUG
        
        

        # MPI stuff
        self.workers = Workers(workers)
                        
        # update 
        self.update()
    
    @property
    def nl(self):
        """Length of the device+contacts"""
        return self._nl

    @nl.setter
    def nl(self, value):
        self._nl = value
        self.update()

    @property
    def nw(self):
        """Width of the device"""
        return self._nw

    @nw.setter
    def nw(self, value):
        self._nw = value
        self.update()

    
             
    ## 
    # Update internally calculated data members 
    # @Parameters:
    #    workers: mpi.world - MPI world communication between worker processes.
    #    
    def update(self):
        # Bounding box
        self.lx = (self.nl-1)*self.a
        self.ly = (self.nw-1)*self.a
        self.xmn =-self.lx/2 - self.a*7/22 + 1/sqrt(2)   # just to make sure that no point is  
        self.xmx = self.lx/2 + self.a*7/22 + 1/sqrt(2)   # at the border
        self.ymn =-self.ly/2 - self.a*7/22 + 1/sqrt(2)
        self.ymx = self.ly/2 + self.a*7/22 + 1/sqrt(2)

   
    def __str__(self):
        msg  = " TI Simulation Parameters:\n"
        msg += "  Device:\n"
        msg += "   Length " + str(self.nl) + " (" + str(self.lx/10.0) + " nm)" 
        msg +=    ", width " + str(self.nw) + " (" + str(self.ly/10.0) + " nm)\n"
        msg += "  k.p parameters:\n"
        msg += "   a = " + str(self.a) + ", dtol = " + str(self.dtol) + ".\n" 
        msg += "   C = " + str(self.C) + ", A2 = " + str(self.A2) + ", K = "
        msg += str(self.K) + ".\n" 
        msg += "  Band Struct parameters:\n"
        msg += "   Dimension = " + str(self.Dim) + ", Fermi Level = " + str(self.mu) + "\n"
        msg += "   Save output at: " + self.OutPath + self.OutFileName + "*.\n" 
        
        return msg
    
    ## 
    # Saves simulation parameters to a file for future use.
    #
    def save(self):
        if (self.workers.IAmMaster()):
            # Create directory if not exist.
            if not os.path.exists(self.OutPath):
                os.makedirs(self.OutPath)  
            
            # Pickle filename
            self.pickleFileName = self.OutFileName + "_simu.pkl"
            output = open(self.OutPath + self.pickleFileName, 'wt')
            pickle.dump(self, output)
            output.close()
            #pass
    
    def __getstate__(self):
        d = dict(self.__dict__)
        del d['workers']
        return d
    
    def __setstate__(self, d):
        self.__dict__.update(d)

## 
# Loads and returns a TiSurfSimuPatams from pickle file
#
def load(fileName):
    pf = open(fileName)
    tisimu = pickle.load(pf)
    pf.close()
    return tisimu


        
