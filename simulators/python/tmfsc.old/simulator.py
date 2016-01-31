"""

"""
from device import Device, nm, AA, EDGE_REFLECT, EDGE_ABSORB

from math import cos, sin, tan, pi, sqrt, atan2
import numpy as np
import numpy.random as rnd

class Simulator(object):
    def __init__(self, dev, mpiworld=None):
        self.dev = dev
        self.max_num_steps = 10000
        self.pts_per_cycles = 100
        self.num_dt_step = 1000
        self.vF = 1E6

        self.dl = 50
        self.nth = 50

        self.mpiworld = mpiworld

    def calc_transmission(self, B, EF, V, cont_num = 0, draw=False, show_progress=True):
        dev = self.dev
        ie = dev.contacts[cont_num]
        if dev.edge_types[ie] != EDGE_ABSORB:
            raise Exception(" Edge # " + str(ie) + " is not a contact.")

        npts = int(dev.edge_lengths[ie]/self.dl)
        inj_pts = dev.get_ptsonedge(cont_num, npts)

        edge_vec = dev.edge_vec[ie]
        th0 = atan2(edge_vec[1], edge_vec[0])

        # reset the bins
        nc = len(dev.contacts)
        for ic in range(nc):
            ie = dev.contacts[ic]
            dev.collected[ie] = 0

        ne = 0
        for ip in range(npts):
            ri = inj_pts[ip, :]*AA

            th = rnd.normal(0, pi/5, self.nth)

            for thi in th:
                if abs(thi) < pi/2.0-pi/20.0:
                    self.calc_trajectory(ri, th0 + thi + pi/2, B, EF, V)
                    ne += 1
                    if draw == True:
                        dev.draw_trajectory(color='b', alpha=0.05)
            if show_progress == True:
                print (str(int(float(ip+1)*100.0/float(npts))) + "%")


        T = np.zeros((nc,nc))
        for ic in range(nc):
            if ic != cont_num:
                ie = dev.contacts[ic]
                T[cont_num][ic] = float(dev.collected[ie])/float(ne)
        return T


    def calc_trajectory(self, ri, thi, B, EF, V):
        r = []
        r.insert(0, ri)
        dev = self.dev
        vF = self.vF
        vi = vF*np.array([cos(thi), sin(thi)])
        wc = vF**2*B/(EF-V)
        dt = abs(2*pi/wc/self.pts_per_cycles)
        dth = wc*dt

        ii = 0
        while ii < self.max_num_steps:

            vf, rf, thf = self.calc_next_state(vi, thi, ri, dth, dt)
            if ii != 0:
                n,x,y = dev.intersects(ri, rf)
                if n != -1:
                    if dev.edge_types[n] == EDGE_REFLECT:
                        dt2 = sqrt((x-ri[0])**2+(y-ri[1])**2)/vF
                        # the following loop might need optimization.
                        while dt2 > 0:
                            dth2 = wc*dt2
                            vf, rf, thf = self.calc_next_state(vi, thi, ri, dth2, dt2)
                            n2,x,y = self.dev.intersects(ri, rf)
                            if n2 == -1:
                                break
                            dt2 = dt2 - dt/self.num_dt_step

                        a1 = sum(vf*dev.edge_vec[n])*dev.edge_vec[n]
                        a2 = vf - a1
                        vf = a1 - a2
                        thf = atan2(vf[1], vf[0])
                    elif dev.edge_types[n] == EDGE_ABSORB:
                        dev.collected[n] = dev.collected[n] + 1
                        r.append(rf)
                        break
            r.append(rf)
            ri = rf
            thi = thf
            vi = vf
            ii += 1
        N = len(r)
        self.dev.trajectory = np.zeros((len(r), 2))
        for ii in range(N):
            self.dev.trajectory[ii,:] = r[ii]/AA

        return rf, thf

    def calc_next_state(self, vi, thi, ri, dth, dt):
        thf = thi + dth
        vf = self.vF*np.array([cos(thf), sin(thf)])
        rf = ri + (vi + vf)/2*dt
        return vf, rf, thf
