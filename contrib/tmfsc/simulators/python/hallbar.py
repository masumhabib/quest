#!/usr/bin/python

from qmicad import greet
from qmicad.tmfsc import Device, Simulator, Trajectory
from qmicad.tmfsc import nm, AA, EDGE_REFLECT, EDGE_ABSORB, EDGE_TRANSMIT
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib import animation

import numpy as np
from math import pi, sqrt, sin, cos, tan
import sys
import os
import shutil
import operator as op
import time

#import boost.mpi as mpi
import mpi

# Constants
vf 			= 1E6		# Fermi velocity

class Bias(object):
    def __init__(self, biasVars = {}):
        self.biasVars = biasVars
        self.sizes = {}
        for key in self.biasVars:
            self.sizes[key] = len(self.biasVars[key]);

    def clear(self):
        self.biasVars = {}
        self.sizes = {}

    def append(self, bias, key):
        self.biasVars[key] = bias
        self.sizes[key] = len(self.biasVars[key]);

    def get(self, ib):
        biases = []
        for key in self.biasVars:
            bias = self.biasVars[key]
            size = self.sizes[key];
            indx = ib % size 
            biases.append(bias[indx])
            ib = ib/size
        return biases

    def __str__(self):
        msg = ""
        for key in self.biasVars:
            msg += "\nBias List " + str(key) + ": "
            msg += '  '.join('%.3f' % bias for bias in self.biasVars[key])
        return msg

    def numBiases(self):
        n = 1
        for key in self.biasVars:
            n = n*self.sizes[key]
        return n

class HallBar(object):
    def __init__(self, mpiworld=None):
        """Constructs a HallBar structure with four contacts."""
        self.fontSize = 20
        self.outDir = './'
        self.outFile = 'trans.npz'
        self.mpiworld = mpiworld
        self.vF = vf
        self.bias = Bias()

        self.clear()
        self.setupDefaults()

    def clear(self):
        """ Clears ourselves """
        # Create device
        self.dev = Device()
        # Create simulator
        self.sim = Simulator(self.dev)
 
    def setupDefaults(self):
        """Sets up a default simulation with a default goemetry."""
        # Fermi level is set to 0, you wont probably need to change this
        self.EF = 0.0	    # Fermi level		

        # Create the goemetry with the following dimensions
        lx = 500.0		# length
        ly = 500.0		# width
        clx = 50.0		# contact width
        cly = 1.0		# contact length
        coffx = lx/10
        dc = (lx - 2.0*coffx - clx) # distance between contacts
        self.dc = dc;
	
        # Create device and add vertices
        self.addPoint(-lx/2, -ly/2)
        self.addPoint(-lx/2+coffx, -ly/2)
        self.addPoint(-lx/2+coffx, -ly/2-cly)
        self.addPoint(-lx/2+coffx+clx, -ly/2-cly)
        self.addPoint(-lx/2+coffx+clx, -ly/2)
        self.addPoint(lx/2-coffx-clx, -ly/2)
        self.addPoint(lx/2-coffx-clx, -ly/2-cly*50)
        self.addPoint(lx/2-coffx, -ly/2-cly*50)
        self.addPoint(lx/2-coffx, -ly/2)
        self.addPoint(lx/2, -ly/2)
        self.addPoint(lx/2, ly/2)
        npts = self.addPoint(-lx/2, ly/2)
        npts = npts + 1
        # add reflecting edges
        for ip in range(npts):
            self.addEdge(ip, (ip+1) % npts);
        # Contacts 1 and 2
        self.setEdgeType(2, EDGE_ABSORB)
        self.setEdgeType(6, EDGE_ABSORB)
        # Contacts 3 and 4
        self.setEdgeType(9, EDGE_ABSORB)
        self.setEdgeType(11, EDGE_ABSORB)

        # Add the transmitting edge
        #ipt1 = self.addPoint(0, -ly/2)
        #ipt2 = self.addPoint(0, ly/2)
        #self.addEdge(ipt1, ipt2, EDGE_TRANSMIT)

        # No gates for default device
        self.gates = [] # gate geometries for drawing

        # default bias setup 
        self.setupBias(1.0, 0.15)

    def addPoint(self, x, y):
        return self.dev.addPoint(np.array([x, y]))

    def addEdge(self, ipt1, ipt2, type = EDGE_REFLECT):
        return self.dev.addEdge(ipt1, ipt2, type);

    def addGate(self, lb, rb, rt, lt, VgRatio = 1):
        gate = (np.array([lb[0], lb[1]]), np.array([rb[0], rb[1]]), 
                np.array([rt[0], rt[1]]), np.array([lt[0], lt[1]]))
        self.gates.append(gate)
        return self.dev.addGate(gate[0], gate[1], gate[2], gate[3]);

    def setEdgeType(self, id, type):
        self.dev.edgeType(id, type)
        
    def setSplitLen(self, len):
        self.dev.setSplitLen(len)
    def getSplitLen(self):
        return self.dev.getSplitLen()

    def setupBias(self, B, V, m = 1, singleResonance = True, 
            Bmax = None, NB = 1, Vmax = None, NV = 1):
        """ Sets up the bias points """
        VV = [] 
        BB = [] 
        if NV == 1:
            VV = np.array([V])
        else:
            assert Vmax is not None, "Vmax is None"
            VV = np.linspace(V, Vmax, NV)

        if NB == 1:   
            BB = np.array([B])
        else:
            assert Bmax is not None, "Vmax is None"
            BB = np.linspace(B, Bmax, NB)
        if singleResonance == True:
            B0 = abs(self.EF-V)/vf/nm/self.dc*2.0
            BB = m*np.array([B0])
        self.bias.clear()
        self.bias.append(BB, 'B')
        self.bias.append(VV, 'V')
 
    def enableDirectCalc(self): 
        self.sim.ParticleType = 1

    def enableCyclotronCalc(self): 
        self.sim.ParticleType = 0

    def numPointsPerCycle(self, npts):
        self.sim.NumPointsPerCycle = npts

    def calcSingleTraject(self, shiftxy=(0,0), shiftth=0, contId=0):
        """ Calculate trajectory for single injection event from the
            midpoint of contId"""
        thi = self.dev.contDirctn(contId) + pi/2 + shiftth
        ri = self.dev.contMidPoint(contId) + 1E-3*self.dev.contNormVect(contId) + np.array([shiftxy[0], shiftxy[1]])
        
        print ("\nInjecting electron from " + str(ri) + " ")
        
        # calculate the trajectory
        biases = self.bias.get(0);
        self.printBias(biases, 0)
        print ""
        if self.dev.NumGates > 0:
            self.trajs = self.sim.calcTraj(ri, thi, biases[0], self.EF, 
                    biases[1:])
        else:
            self.trajs = self.sim.calcTraj(ri, thi, biases[0], self.EF, 
                    biases[1])
    
    def calcSingleTrans(self, dl=5, nth=50, saveTrajectory=False, 
            contId = 0):
        """ Transmission for single B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        biases = self.bias.get(0);
        self.printBias(biases, 0)
        if self.dev.NumGates > 0:
            self.T,self.trajs = self.sim.calcTrans(biases[0], self.EF, 
                    biases[1:], contId, saveTrajectory)
        else:
            self.T,self.trajs = self.sim.calcTrans(biases[0], self.EF, 
                    biases[1], contId, saveTrajectory)
            self.printTrans(contId, self.T)
        
        return self.T
    
    def calcAllTrans(self, dl=5, nth=50, contId=0):
        """ Transmission for all B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        
        nconts = self.dev.numConts()
        npts = self.bias.numBiases()
        self.T = np.zeros((npts, nconts, nconts))
        
        # MPI stuff
        ncpu = self.mpiworld.size
        icpu = self.mpiworld.rank
        quo = npts/ncpu
        rem = npts%ncpu
        my_start = icpu*quo
        if icpu < rem:
            my_start += icpu
        else:
            my_start += rem
        my_end = my_start + quo
        if icpu < rem:
            my_end += 1

        print "\nTransmission loop ..."
        for ipt in range(my_start, my_end):
            biases = self.bias.get(ipt)
            self.printBias(biases, ipt)
            if self.dev.NumGates > 0:
                T,self.trajs = self.sim.calcTrans(biases[0], self.EF, 
                        biases[1:], False, contId)
            else:
                T,self.trajs = self.sim.calcTrans(biases[0], self.EF, 
                        biases[1], False, contId)
                self.printTrans(contId, T)
            self.T[ipt,:,:] = T
 
        if self.mpiworld.rank == 0:
            self.T = mpi.reduce(self.mpiworld, self.T, op.add, 0) 
        else:
            mpi.reduce(self.mpiworld, self.T, op.add, 0) 
 
    def drawGeom(self):
        """ Draws the device outline """
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(111)
       
        nedges = self.dev.numEdges()
        edgeVec = self.dev.edgeVect(0)
        pt1 = self.dev.edgeMidPoint(0) - edgeVec/2.0
        for ie in range(nedges):                                    
            if self.dev.edgeType(ie) == EDGE_ABSORB:
                width = 4.0
            else:
                width = 1.5 
            edgeVec = self.dev.edgeVect(ie)
            pt2 = pt1 + edgeVec
            X = np.array([pt1[0], pt2[0]])
            Y = np.array([pt1[1], pt2[1]])
            self.axes.plot(X, Y, 'r-', linewidth=width) 
            pt1 = pt2
        # draw gates
        codes = [Path.MOVETO,
            Path.LINETO,
            Path.LINETO,
            Path.LINETO,
            Path.CLOSEPOLY,
         ]
        for gate in self.gates:
            gate = gate + (np.array([0,0]),)
            path = Path(gate, codes)
            patch = patches.PathPatch(path, alpha=0.5, facecolor='yellow', 
                lw=1)
            self.axes.add_patch(patch)

        self.axes.set_aspect('equal', 'datalim')
        self.axes.set_xlabel('x (nm)')
        self.axes.set_ylabel('y (nm)')

        ncnts = self.dev.numConts()
        for ic in range(ncnts):
            pt = self.dev.contMidPoint(ic) - 30*self.dev.contNormVect(ic)
            self.axes.text(pt[0], pt[1], str(ic+1), fontsize=self.fontSize)


    def drawTrajectory(self, color=None, alpha=1.0, width=2.0):
        for traj in self.trajs:
            alpha = traj.occupation
            if color is None:
                if alpha > 1:
                    alpha = 1
                self.axes.plot(traj.path[:, 0], traj.path[:, 1],
                    linewidth=width, alpha=alpha)
            else:
                 self.axes.plot(traj.path[:, 0], self.traj.path[:, 1], 
                         linewidth=width, color=color, alpha=alpha)
 
    def animate(self, filename=None):
	    self._start_animation()
        #if filename is not None:
	    #    self.dev.save_animation(filename)

    def showPlot(self):
        plt.show()


    def saveTraj(self, file_name):
        plt.savefig(file_name, dpi=300)

    def save(self):
        """ Saves results """
        if self.mpiworld.rank == 0:
            if not os.path.exists(self.outDir):
                os.makedirs(self.outDir)
            np.savez_compressed(self.outDir+self.outFile, T=self.T, B=self.B, V=self.V) 


    def printBiasList(self):
        print (self.bias)

    def printTraj(self):
        print("")
        print("Trajectory:")
        print(self.trajs)
    
    def printBias(self, biases, ipt):
        msg = "Bias# " + str(ipt) + ": "
        biasIndx = range(1, len(biases));
        msg += "B={0:.3f}".format(biases[0])
        msg += '  '.join('  V{0}={1:.3f}'.format(ib, biases[ib]) for ib in biasIndx)
        print msg,
 
    def printTrans(self, contId, T):
        msg = " ==>"
        for ic in range(self.dev.numConts()):
            if ic != contId:
                msg += '  T{0}{1}={2:.3f}'.format(contId, ic, T[contId][ic])
        print (msg)
 
    def banner(self):
        print(greet())
        print("\n ***  Running semiclassical analysis for graphene ... ")

#    def _getGateVoltages(self, V):
#        VGs = []
#        for VgRatio in self.VgRatios:
#            VGs.append(self.VGG + V*VgRatio)
#        return VGs

    def _start_animation(self): 
        fig = self.fig                                                          
        ax = self.axes                                                          
        self.line, = ax.plot([], [], 'ro')                                      
        self.ani = animation.FuncAnimation(fig, self._set_anim_pos, 
            self._anim_data, blit=False, interval=10, repeat=True) 
        plt.show()    

    def _anim_data(self):
        traj = self.trajs[0]
        for xy in traj:
            yield xy[0], xy[1]

    def _set_anim_pos(self, data): 
        x, y = data[0], data[1]                                                 
        self.line.set_data(x, y)                                                
        return self.line      

"""
The main() function.
"""
def main(argv = None):
    if argv is None:
        argv = sys.argv

    if len(argv) > 1:
        start = time.time()

        simu_file = argv[1];
        if not os.path.exists(simu_file):
            raise Exception("File " + simu_file + " not found")

        hallbar = HallBar(mpi.world)
        hallbar.banner()
        exec(open(simu_file).read(), globals(), locals())

        if not os.path.exists(hallbar.outDir):
            os.makedirs(hallbar.outDir)
        shutil.copy2(simu_file, hallbar.outDir+simu_file)
        
        elapsed = time.time() - start
        if mpi.world.rank == 0:
            print("\n ***  RUNTIME: " + str(int(elapsed)) + " s.")

    else:
        usage  = " Usage: hallbar.py simu.py\n"
        usage += "   simu.py: The python file defining the simulation parameters.\n"
        print(usage)

    return 0


"""
Entry point.
"""
if __name__ == "__main__":
    sys.exit(main())






