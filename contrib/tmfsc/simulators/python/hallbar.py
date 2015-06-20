#!/usr/bin/python

from qmicad.tmfsc import Device, Simulator
from qmicad.tmfsc import nm, AA, EDGE_REFLECT, EDGE_ABSORB, EDGE_TRANSMIT
import matplotlib.pyplot as plt
from matplotlib import animation

import numpy as np
from math import pi
import sys
import os
import shutil
import operator as op
import time

#import boost.mpi as mpi
import mpi

# Constants
TWO_CONTS  	= 2 		# two contacts
FOUR_CONTS 	= 4   		# four contacts
vf 			= 1E6		# Fermi velocity

class HallBar(object):
    def __init__(self, mpiworld=None):
        # simulation setup
        self.num_contacts = FOUR_CONTS

        # bias
        self.EF = 0.0	    # Fermi level		
        self.Bmin = 1.0	    # Magnetic field
        self.Vmin = 0.15	# Gate bias
        self.NB = 1         # size of self.B
        self.NV = 1         # size of self.V
        self.m  = 0.99      # Resonance number (optional)

        # dimensions
        self.lx = 500.0		# length
        self.ly = 500.0		# width
        self.clx = 50.0		# contact width
        self.cly = 1.0		    # contact length

        self.fontsize = 20

        self.out_dir = './'
        self.out_file = 'trans.npz'

        self.mpiworld = mpiworld

    def setup_geom(self):
        """ Creates the goemetry """
        lx = self.lx
        ly = self.ly
        clx = self.clx
        cly = self.cly
        self.coffx = self.lx/10
        self.dc = (self.lx - 2.0*self.coffx - self.clx) # distance between contacts
	
        coffx = self.coffx
        dc = self.dc

        # Create device and add vertices
        self.dev = Device()
        self.dev.addPoint(np.array([-lx/2, -ly/2]))
        self.dev.addPoint(np.array([-lx/2+coffx, -ly/2]))
        self.dev.addPoint(np.array([-lx/2+coffx, -ly/2-cly]))
        self.dev.addPoint(np.array([-lx/2+coffx+clx, -ly/2-cly]))
        self.dev.addPoint(np.array([-lx/2+coffx+clx, -ly/2]))
        self.dev.addPoint(np.array([lx/2-coffx-clx, -ly/2]))
        self.dev.addPoint(np.array([lx/2-coffx-clx, -ly/2-cly*50]))
        self.dev.addPoint(np.array([lx/2-coffx, -ly/2-cly*50]))
        self.dev.addPoint(np.array([lx/2-coffx, -ly/2]))
        self.dev.addPoint(np.array([lx/2, -ly/2]))
        self.dev.addPoint(np.array([lx/2, ly/2]))
        self.dev.addPoint(np.array([-lx/2, ly/2]))
        self.dev.addPoint(np.array([-lx/2, -ly/2]))

        # Contacts 1 and 2
        self.dev.edgeType(2, EDGE_ABSORB)
        self.dev.edgeType(6, EDGE_ABSORB)

        # Csontacts 3 and 4
        if self.num_contacts == FOUR_CONTS:
            self.dev.edgeType(9, EDGE_ABSORB)
            self.dev.edgeType(11, EDGE_ABSORB)

        # Create simulator
        self.sim = Simulator(self.dev)
        self.vF = vf

    def setup_bias(self):
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


    def calc_single_traject(self, shiftxy=(0,0), shiftth=0):
        """ Calculate trajectory for single injection event """
        thi = pi/2 + shiftth
        ri = np.array([-self.lx/2+self.coffx+self.clx/2+shiftxy[0]+1E-3, -self.ly/2-self.cly+shiftxy[1]+1E-3])
        
        print ("Injecting electron from " + str(ri) + " ")
        
        # calculate the trajectory
        self.trajs = [];
        self.trajs.append(self.sim.calcTraj(ri, thi, self.B[0], self.EF, 
            self.V[0]))
    
    def calc_single_trans(self, dl=50, nth=50, saveTrajectory=False, 
            contId = 0):
        """ Transmission for single B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        T,self.trajs = self.sim.calcTrans(self.B[0], self.EF, self.V[0], 
                contId, saveTrajectory)
        
        T12 = ' T12 = ' + '{0:.2f}'.format(T[0][1])
        if self.num_contacts == FOUR_CONTS:
	        T13 = ' T13 = ' + '{0:.2f}'.format(T[0][2])
	        T14 = ' T14 = ' + '{0:.2f}'.format(T[0][3])
        else:
	        T13 = ''
	        T14 = ''
        print ('B = {0:.3f}'.format(self.B[0]) + ' V = {0:.3f}'.format(self.V[0]) + T12+T13+T14)
        
        return T
    
    def calc_all_trans(self, dl=50, nth=50):
        """ Transmission for all B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        
        self.T12 = np.zeros((self.NB, self.NV))
        self.T13 = np.zeros((self.NB, self.NV))
        self.T14 = np.zeros((self.NB, self.NV))
        
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

        for ipt in range(my_start, my_end):
            ib = ipt/self.NV
            iv = ipt%self.NV
            T = self.sim.calc_transmission(self.B[ib], self.EF, self.V[iv], draw=False, show_progress=False)
            
       
            self.T12[ib][iv] = T[0][1]
            T12 = ' T12 = ' + '{0:.2f}'.format(T[0][1])

            if self.num_contacts == FOUR_CONTS:
                self.T13[ib][iv] = T[0][2]
                self.T14[ib][iv] = T[0][3]
 
                T13 = ' T13 = ' + '{0:.2f}'.format(T[0][2])
                T14 = ' T14 = ' + '{0:.2f}'.format(T[0][3])
            else:
                T13 = ''
                T14 = ''
            print ('B = {0:.3f}'.format(self.B[ib]) + ' V = {0:.3f}'.format(self.V[iv]) + T12+T13+T14)
 
        if self.mpiworld.rank == 0:
            self.T12 = mpi.reduce(self.mpiworld, self.T12, op.add, 0) 
            self.T13 = mpi.reduce(self.mpiworld, self.T13, op.add, 0) 
            self.T14 = mpi.reduce(self.mpiworld, self.T14, op.add, 0) 
        else:
            mpi.reduce(self.mpiworld, self.T12, op.add, 0) 
            mpi.reduce(self.mpiworld, self.T13, op.add, 0) 
            mpi.reduce(self.mpiworld, self.T14, op.add, 0) 
 
    def draw_geom(self):
        """ Draws the device outline """
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(111)
       
        nedges = self.dev.numEdges()
        pt1 = np.array([-self.lx/2, -self.ly/2])
        for ie in range(0, nedges):                                    
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

        xc1 = (-self.lx/2+self.coffx+self.clx/2-10)
        yc1 = (-self.ly/2-self.cly-60)
        self.axes.text(xc1, yc1, '1', fontsize=self.fontsize)
        xc2 = (self.lx/2-self.coffx-self.clx+10)
        yc2 = (-self.ly/2-self.cly*50-60)
        self.axes.text(xc2, yc2, '2', fontsize=self.fontsize)
        if self.num_contacts == FOUR_CONTS:
            xc3 = (self.lx/2+30.0)
            yc3 = (-10.0)
            self.axes.text(xc3, yc3, '3', fontsize=self.fontsize)
            xc4 = (-self.lx/2-50.0)
            yc4 = (-10.0)
            self.axes.text(xc4, yc4, '4', fontsize=self.fontsize)


    def draw_trajectory(self, color=None, alpha=1.0, width=2.0):
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

    def show_plot(self):                                                        
        plt.show()


    def save_trajectory(self, file_name):
        plt.savefig(file_name, dpi=300)

    def save(self):
        """ Saves results """
        if self.mpiworld.rank == 0:
            if not os.path.exists(self.out_dir):
                os.makedirs(self.out_dir)
            np.savez_compressed(self.out_dir+self.out_file, T12=self.T12, T13=self.T13, T14=self.T14, B=self.B, V=self.V) 


    def print_bias(self):
        print("")
        print("Bias List:")
        print(self.V)
        print("")
        print("B Field List:")
        print(self.B)

    def print_traject(self):
        print("")
        print("Trajectory:")
        print(self.trajectory)

    def banner(self):
        pass

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
        exec(open(simu_file).read(), globals(), locals())

        if not os.path.exists(hallbar.out_dir):
            os.makedirs(hallbar.out_dir)
        shutil.copy2(simu_file, hallbar.out_dir+simu_file)
        
        elapsed = time.time() - start
        if mpi.world.rank == 0:
            print(' RUNTIME: ' + str(int(elapsed)) + ' sec')

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






