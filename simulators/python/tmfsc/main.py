#!/opt/local/bin/python

import device
from device import Device, nm, AA
from simulator import Simulator
from intersect import intersects, intersection, test
import matplotlib.pyplot as plt

import numpy as np
from math import pi

B = 0.25
EF = 0.0
V = 0.15

w = 50000
l = 100000
lc = l/100
wc = w/100
lcoff = l/10

dev = Device()

dev.add_pts(np.array([[-l/2, -w/2],
                      [-l/2+lcoff, -w/2],
                      [-l/2+lcoff, -w/2-wc],
                      [-l/2+lcoff+lc, -w/2-wc],
                      [-l/2+lcoff+lc, -w/2],
                      [l/2-lcoff-lc, -w/2],
                      [l/2-lcoff-lc, -w/2-wc],
                      [l/2-lcoff, -w/2-wc],
                      [l/2-lcoff, -w/2],
                      [l/2, -w/2],
                      [l/2, w/2],
                      [-l/2, w/2],
                      [-l/2, -w/2]
                      ]))

dev.set_edge_type(2, device.EDGE_ABSORB)
dev.set_edge_type(6, device.EDGE_ABSORB)

# dev.set_edge_type(9, device.EDGE_ABSORB)
# dev.set_edge_type(11, device.EDGE_ABSORB)

sim = Simulator(dev)
# thi = pi/2;
# ri = [0, -w/2]*nm;
cont1 = np.array([-l/2+lcoff+lc/2, -w/2-wc])*AA
thi = pi/2
ri = cont1

[r, th] = sim.calc_trajectory(ri, thi, B, EF, V)
# # x = sim.trajectory[:, 0]
# # y = sim.trajectory[:, 1]
#
dev.draw_geom()
dev.draw_trajectory()
plt.savefig('reflecting1.png', dpi=300)
plt.show()
# dev.start_animation()
# dev.save_animation('reflection.mp4')

# l1 = dev.edges[4]/10
# l2 = np.array([[-l/2+lcoff+lc/2, -w/2+w/4], [-l/2+lcoff+lc/2, -w/2-w/4]])/10
# test(l1, l2)
