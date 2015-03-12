#!/usr/bin/python

import device
from device import Device, nm, AA
from simulator import Simulator
from intersect import intersects, intersection, test
import matplotlib.pyplot as plt

import numpy as np
from math import pi
import sys
import os
import shutil

import boost.mpi as mpi

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
        self.Bmin = 1	    # Magnetic field
        self.Vmin = 0.15	# Gate bias
        self.NB = 1         # size of self.B
        self.NV = 1         # size of self.V
        self.m  = 0.99      # Resonance number (optional)

        # dimensions
        self.lx = 5000		# length
        self.ly = 5000		# width
        self.clx = 500		# contact width
        self.cly = 10		# contact length

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
        self.dc = (self.lx - 2.0*self.coffx - self.clx)*AA # distance between contacts
	
        coffx = self.coffx
        dc = self.dc

        # Create device and add vertices
        self.dev = Device()
        self.dev.add_pts(np.array([[-lx/2, -ly/2],
				              	   [-lx/2+coffx, -ly/2],
				              	   [-lx/2+coffx, -ly/2-cly],
				              	   [-lx/2+coffx+clx, -ly/2-cly],
				                   [-lx/2+coffx+clx, -ly/2],
				                   [lx/2-coffx-clx, -ly/2],
				                   [lx/2-coffx-clx, -ly/2-cly*50],
				                   [lx/2-coffx, -ly/2-cly*50],
				                   [lx/2-coffx, -ly/2],
				                   [lx/2, -ly/2],
				                   [lx/2, ly/2],
				                   [-lx/2, ly/2],
				                   [-lx/2, -ly/2]
				                  ]))

        # Contacts 1 and 2
        self.dev.set_edge_type(2, device.EDGE_ABSORB)
        self.dev.set_edge_type(6, device.EDGE_ABSORB)

        # Csontacts 3 and 4
        if self.num_contacts == FOUR_CONTS:
            self.dev.set_edge_type(9, device.EDGE_ABSORB)
            self.dev.set_edge_type(11, device.EDGE_ABSORB)

        # Create simulator
        self.sim = Simulator(self.dev, self.mpiworld)
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
                B0 = abs(self.EF-self.V)/vf/self.dc*2.0
                self.B = self.m*B0
        else:
            self.B = np.linspace(self.Bmin, self.Bmax, self.NB)

    def draw_geom(self):
        """ Draws the device outline """
        self.dev.draw_geom()
        xc1 = (-self.lx/2+self.coffx+self.clx/2-100)/10
        yc1 = (-self.ly/2-self.cly-600)/10
        self.dev.axes.text(xc1, yc1, '1', fontsize=self.fontsize)
        xc2 = (self.lx/2-self.coffx-self.clx+100)/10
        yc2 = (-self.ly/2-self.cly*50-600)/10
        self.dev.axes.text(xc2, yc2, '2', fontsize=self.fontsize)
        if self.num_contacts == FOUR_CONTS:
            xc3 = (self.lx/2+300.0)/10
            yc3 = (-100.0)/10
            self.dev.axes.text(xc3, yc3, '3', fontsize=self.fontsize)
            xc4 = (-self.lx/2-500.0)/10
            yc4 = (-100.0)/10
            self.dev.axes.text(xc4, yc4, '4', fontsize=self.fontsize)


    def calc_single_traject(self, filename=None, animate=False, shiftxy=(0,0), shiftth=0):
        """ Calculate trajectory for single injection event """
        thi = pi/2 + shiftth
        ri = np.array([-self.lx/2+self.coffx+self.clx/2+shiftxy[0], -self.ly/2-self.cly+shiftxy[1]])*AA
        
        
        [r, th] = self.sim.calc_trajectory(ri, thi, self.B[0], self.EF, self.V[0])

        self.dev.draw_trajectory()

        if animate == True:
	        self.dev.start_animation()

        self.dev.show_plot()
        		
        if filename is not None:
            if animate == True:
	            self.dev.save_animation(filename)
            else:
                self.dev.save_trajectory(filename)


    def calc_single_trans(self, dl=50, nth=50, draw=False):
        """ Transmission for single B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        T = self.sim.calc_transmission(self.B[0], self.EF, self.V[0], draw=draw)
        
        T12 = ' T12 = ' + '{0:.2f}'.format(T[0][1])
        if self.num_contacts == FOUR_CONTS:
	        T13 = ' T13 = ' + '{0:.2f}'.format(T[0][2])
	        T14 = ' T14 = ' + '{0:.2f}'.format(T[0][3])
        else:
	        T13 = ''
	        T14 = ''
        print ('B = {0:.3f}'.format(self.B[0]) + ' V = {0:.3f}'.format(self.V[0]) + T12+T13+T14)
        
        return T
    
    def calc_all_traject(self, filename=None, dl=50, nth=50):  
        """ Calculates all the trajectories from one contact """      
        T = self.calc_single_trans(dl=dl, nth=nth, draw=True)

        T12 = '$T_{12}= ' + '{0:.2f}$'.format(T[0][1])
        if self.num_contacts == FOUR_CONTS:
	        T13 = '$T_{13}= ' + '{0:.2f}$'.format(T[0][2])
	        T14 = '$T_{14}= ' + '{0:.2f}$'.format(T[0][3])
        else:
	        T13 = ""
	        T14 = ""

        xt = -3500/10
        yt = (self.ly/2+200)/10
        self.dev.axes.text(xt, yt, T12 + " " + T13 + " " + T14, fontsize=self.fontsize)

        if filename is not None:
	        self.dev.save_trajectory(filename)

        self.dev.show_plot()

    def calc_all_trans(self, dl=50, nth=50):
        """ Transmission for all B and V """
        self.sim.dl = dl
        self.sim.nth = nth
        
        self.T12 = np.zeros((self.NB, self.NV))
        self.T13 = np.zeros((self.NB, self.NV))
        self.T14 = np.zeros((self.NB, self.NV))

        ib = 0
        for B in self.B:
            iv = 0
            for V in self.V:
                T = self.sim.calc_transmission(B, self.EF, V, draw=False, show_progress=False)
                self.T12[ib][iv] = T[0][1]
                self.T13[ib][iv] = T[0][2]
                self.T14[ib][iv] = T[0][3]
        
                # User feedback
                T12 = ' T12 = ' + '{0:.2f}'.format(T[0][1])
                if self.num_contacts == FOUR_CONTS:
	                T13 = ' T13 = ' + '{0:.2f}'.format(T[0][2])
	                T14 = ' T14 = ' + '{0:.2f}'.format(T[0][3])
                else:
	                T13 = ''
	                T14 = ''
                print ('B = {0:.3f}'.format(self.B[ib]) + ' V = {0:.3f}'.format(self.V[iv]) + T12+T13+T14)
 
                iv += 1
            ib += 1

    def save(self):
        """ Saves results """
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        np.savez_compressed(self.out_dir+self.out_file, T12=self.T12, T13=self.T13, T14=self.T14, B=self.B, V=self.V) 





"""
The main() function.
"""
def main(argv = None):
    if argv is None:
        argv = sys.argv

    if len(argv) > 1:
        simu_file = argv[1];
        if not os.path.exists(simu_file):
            raise Exception("File " + simu_file + " not found")

        hallbar = HallBar(mpi.world)
        exec(open(simu_file).read(), globals(), locals())

        if not os.path.exists(hallbar.out_dir):
            os.makedirs(hallbar.out_dir)
        shutil.copy2(simu_file, hallbar.out_dir+simu_file)
        

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






