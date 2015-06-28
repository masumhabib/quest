#!/usr/bin/python

from qmicad import greet
from qmicad.tmfsc import Device, Simulator
from qmicad.tmfsc import nm, AA, EDGE_REFLECT, EDGE_ABSORB, EDGE_TRANSMIT
import matplotlib.pyplot as plt
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

class HallBar(object):
    def __init__(self, mpiworld=None):
        """Constructs a HallBar structure with four contacts."""

        self.fontSize = 20

        self.outDir = './'
        self.outFile = 'trans.npz'

        self.mpiworld = mpiworld
       
        self.vF = vf

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
        # bias
        self.EF = 0.0	    # Fermi level		
        self.Bmin = 1.0	    # Magnetic field
        self.Vmin = 0.15	# Gate bias
        self.NB = 1         # size of self.B
        self.NV = 1         # size of self.V
        self.m  = 0.99      # Resonance number (optional)

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
        self.addPoint(-lx/2, ly/2)
        self.addPoint(-lx/2, -ly/2)

        # Contacts 1 and 2
        self.setEdgeType(2, EDGE_ABSORB)
        self.setEdgeType(6, EDGE_ABSORB)
        # Csontacts 3 and 4
        self.setEdgeType(9, EDGE_ABSORB)
        self.setEdgeType(11, EDGE_ABSORB)

    def addPoint(self, x, y):
        self.dev.addPoint(np.array([x, y]))

    def setEdgeType(self, id, type):
        self.dev.edgeType(id, type)

    def setupBias(self):
        """ Sets the bias points """
        if self.NV == 1:
            self.V = np.array([self.Vmin])
        else:
            self.V = np.linspace(self.Vmin, self.Vmax, self.NV)

        if self.NB == 1:        
            if self.m <= 0:
                self.B = np.array([self.Bmin])
            else:
                B0 = abs(self.EF-self.V)/vf/nm/self.dc*2.0
                self.B = self.m*B0
        else:
            self.B = np.linspace(self.Bmin, self.Bmax, self.NB)


    def calcSingleTraject(self, shiftxy=(0,0), shiftth=0, contId=0):
        """ Calculate trajectory for single injection event from the
            midpoint of contId"""
        thi = self.dev.contDirctn(contId) + pi/2 + shiftth
        ri = self.dev.contMidPoint(contId) + 1E-3*self.dev.contNormVect(contId) + np.array([shiftxy[0], shiftxy[1]])
        
        print ("Injecting electron from " + str(ri) + " ")
        
        # calculate the trajectory
        self.trajs = [];
        self.trajs.append(self.sim.calcTraj(ri, thi, self.B[0], self.EF, 
            self.V[0]))
    
    def calcSingleTrans(self, dl=5, nth=50, saveTrajectory=False, 
            contId = 0):
        """ Transmission for single B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        self.T,self.trajs = self.sim.calcTrans(self.B[0], self.EF, self.V[0], 
                contId, saveTrajectory)
        self.printTrans(self.T, contId, self.B[0], self.V[0])
        
        return self.T
    
    def calcAllTrans(self, dl=50, nth=50, contId=0):
        """ Transmission for all B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        
        nconts = self.dev.numConts()
        self.T = np.zeros((self.NB, self.NV, nconts, nconts))
        
        # MPI stuff
        npts = self.NB*self.NV
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

        print("  B         V        EF        T        T         T  ")
        for ipt in range(my_start, my_end):
            ib = ipt/self.NV
            iv = ipt%self.NV
            #print("B = {0:.3f} V = {1:.3f} EF = {2:.3f}".format(self.B[ib], self.V[iv], self.EF))
            T,self.trajs = self.sim.calcTrans(self.B[ib], self.EF, self.V[iv], 
                False, contId)
            self.printTrans(T, contId, self.B[ib], self.V[iv])
            self.T[ib,iv,:,:] = T
 
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

        self.axes.set_aspect('equal', 'datalim')                                
        self.axes.set_xlabel('x (nm)')                                          
        self.axes.set_ylabel('y (nm)')       

        ncnts = self.dev.numConts();
        for ic in range(ncnts):
            pt = self.dev.contMidPoint(ic) - 30*self.dev.contNormVect(ic)
            self.axes.text(pt[0], pt[1], str(ic+1), fontsize=self.fontSize)


    def drawTrajectory(self, color=None, alpha=1.0, width=2.0):
        for traj in self.trajs:
            if color is None:
                self.axes.plot(traj[:, 0], traj[:, 1],
                        linewidth=width, alpha=alpha)
            else:
                self.axes.plot(traj[:, 0], self.traj[:, 1], 
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


    def printBias(self):
        print("\nBias List:")
        print(str(self.V) + "\n")
        print("B Field List:")
        print(str(self.B) + "\n")

    def printTraj(self):
        print("")
        print("Trajectory:")
        print(self.trajs)
    
    def printTrans(self, T, contId, B=0, V=0):
        msg  = '{0:.3f}'.format(B)
        msg += '     {0:.3f}'.format(V)
        for ic in range(self.dev.numConts()):
            if ic != contId:
                #msg += ' ' + str(contId+1) + str(ic+1)
                msg += '    {0:.3f}'.format(T[contId][ic])
        print (msg)
 
    def banner(self):
        print(greet())
        print("\n ***  Running semiclassical analysis for graphene ... ")

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






