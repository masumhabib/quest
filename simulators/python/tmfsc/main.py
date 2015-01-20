#!/opt/local/bin/python

import device
from device import Device, nm, AA
from simulator import Simulator
from intersect import intersects, intersection, test
import matplotlib.pyplot as plt

import numpy as np
from math import pi

four_cont = True
four_cont = False

# bias
n = 1.5
EF = 0.0
V = 0.15

w = 5000
l = 10000
lc = 500
wc = 10
lcoff = l/10

dc = (l - 2.0*lcoff - lc)*AA
vf = 1E6


# resonance condition for B
B0 = abs(EF-V)/vf/dc*2.0
B = n*B0

dev = Device()

dev.add_pts(np.array([[-l/2, -w/2],
                      [-l/2+lcoff, -w/2],
                      [-l/2+lcoff, -w/2-wc],
                      [-l/2+lcoff+lc, -w/2-wc],
                      [-l/2+lcoff+lc, -w/2],
                      [l/2-lcoff-lc, -w/2],
                      [l/2-lcoff-lc, -w/2-wc*50],
                      [l/2-lcoff, -w/2-wc*50],
                      [l/2-lcoff, -w/2],
                      [l/2, -w/2],
                      [l/2, w/2],
                      [-l/2, w/2],
                      [-l/2, -w/2]
                      ]))

dev.set_edge_type(2, device.EDGE_ABSORB)
dev.set_edge_type(6, device.EDGE_ABSORB)

if four_cont == True:
    dev.set_edge_type(9, device.EDGE_ABSORB)
    dev.set_edge_type(11, device.EDGE_ABSORB)

sim = Simulator(dev)

dev.draw_geom()
fontsize = 20
xc1 = (-l/2+lcoff+lc/2-100)/10
yc1 = (-w/2-wc-600)/10
dev.axes.text(xc1, yc1, '1', fontsize=fontsize)
xc2 = (l/2-lcoff-lc+100)/10
yc2 = (-w/2-wc*50-600)/10
dev.axes.text(xc2, yc2, '2', fontsize=fontsize)
if four_cont == True:
    xc3 = (l/2+300.0)/10
    yc3 = (-100.0)/10
    dev.axes.text(xc3, yc3, '3', fontsize=fontsize)
    xc4 = (-l/2-500.0)/10
    yc4 = (-100.0)/10
    dev.axes.text(xc4, yc4, '4', fontsize=fontsize)


# # # Single trajectory
# thi = pi/2
# cont1 = np.array([-l/2+lcoff+lc/2, -w/2-wc])*AA
# thi = pi/2
# ri = cont1
# [r, th] = sim.calc_trajectory(ri, thi, B, EF, V)

# dev.draw_trajectory()


# Transmission
sim.dl = 50
sim.nth = 50
T = sim.calc_transmission(B, EF, V, draw=True)
print T
xt = -3500/10
yt = (w/2+200)/10
T12 = '$T_{12}= ' + '{0:.2f}$'.format(T[0][1])
if four_cont == True:
    T13 = '$T_{13}= ' + '{0:.2f}$'.format(T[0][2])
    T14 = '$T_{14}= ' + '{0:.2f}$'.format(T[0][3])
else:
    T13 = ""
    T14 = ""
dev.axes.text(xt, yt, T12 + " " + T13 + " " + T14, fontsize=fontsize)

if four_cont == True:
    dev.save_trajectory('four_cont1.png')
else:
    dev.save_trajectory('two_cont1.png')

dev.show_plot()


# plt.savefig('reflecting1.png', dpi=300)
# plt.show()
# dev.start_animation()
# dev.save_animation('reflection.mp4')

# l1 = dev.edges[4]/10
# l2 = np.array([[-l/2+lcoff+lc/2, -w/2+w/4], [-l/2+lcoff+lc/2, -w/2-w/4]])/10
# test(l1, l2)

# pts = dev.get_ptsonedge(0, 10)
# plt.plot(pts[:,0]/10, pts[:,1]/10, 'o')
# plt.show()








