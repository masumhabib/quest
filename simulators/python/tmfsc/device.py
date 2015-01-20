
# from shapely.geometry import Point, LineString
from intersect import intersects, intersection

from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

nm = 1E-9   # nanometer
AA = 1E-10   # Angstrom
EDGE_REFLECT = 1
EDGE_ABSORB = 2
EDGE_TRANSMIT = 3


class Device(object):
    def __init__(self):
        self.pts = {}
        self.edges = {}
        self.edge_vec = {}
        self.edge_types = {}
        self.edge_lengths = {}
        self.contacts = {}
        self.collected = {}
        self.trajectory = np.array([])

    # dimensions are in angstroms
    def add_pts(self, pts):
        n = len(self.pts)
        for i in range(len(pts)):
            ip = i + n
            pt = pts[i]
            self.pts[ip] = pt

            if ip > 0:
                pt_im1 = self.pts[i-1]
                self.edges[ip-1] = np.array([pt_im1, pt])
                dx = pt[0] - pt_im1[0]
                dy = pt[1] - pt_im1[1]
                r = sqrt(dx*dx+dy*dy)
                self.edge_lengths[ip-1] = r
                self.edge_vec[ip-1] = np.array([dx, dy])/r
                self.edge_types[ip-1] = EDGE_REFLECT

    def set_edge_type(self, indx, edge_type):
        self.edge_types[indx] = edge_type
        self.collected[indx] = 0
        if edge_type == EDGE_ABSORB:
            self.contacts[len(self.contacts)] = indx

    def intersects(self, pti, ptf):
        pti = pti/AA
        ptf = ptf/AA
        for i in range(len(self.edges)):
            edge = self.edges[i]
            if intersects(edge[0, :], edge[1,:], pti, ptf, False):
                pt = intersection(edge[0, :], edge[1,:], pti, ptf)
                if pt:
                    return i, pt[0]*AA, pt[1]*AA
        return -1, 0, 0

    # Using shapely
    # def intersects(self, pti, ptf):
    #     pti = pti/AA
    #     ptf = ptf/AA
    #     line = LineString([(pti[0], pti[1]), (ptf[0], ptf[1])])
    #
    #     i = 0
    #     for edge in self.edges:
    #         line2 = LineString([(edge[0, 0], edge[0, 1]), (edge[1, 0], edge[1, 1])])
    #         if line.intersects(line2):
    #             pt = line.intersection(line2)
    #             return i, pt.x*AA, pt.y*AA
    #         i += 1
    #     return -1, 0, 0

    def get_ptsonedge(self, cont_num, npts):
        edge = self.edges[self.contacts[cont_num]]
        pt1 = edge[0]
        pt2 = edge[1]
        pts = np.zeros((npts, 2))
        x = np.linspace(pt1[0], pt2[0], npts+2)
        y = np.linspace(pt1[1], pt2[1], npts+2)

        pts[:, 0] = pts[:, 0] + x[1:len(x)-1]
        pts[:, 1] = pts[:, 1] + y[1:len(y)-1]

        return pts

    def draw_geom(self):
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(111)
        for ie in range(0, len(self.edges)):
            edge = self.edges[ie]
            if self.edge_types[ie] == EDGE_ABSORB:
                width = 4.0
            else:
                width = 1.5
            self.axes.plot(edge[:, 0]/10, edge[:, 1]/10, 'r-', linewidth=width)
        self.axes.set_aspect('equal', 'datalim')
        self.axes.set_xlabel('x (nm)')
        self.axes.set_ylabel('y (nm)')

    def draw_trajectory(self, color=None, alpha=1.0, width=2.0):
        if color is None:
            self.axes.plot(self.trajectory[:, 0]/10, self.trajectory[:, 1]/10, linewidth=width, alpha=alpha)
        else:
            self.axes.plot(self.trajectory[:, 0]/10, self.trajectory[:, 1]/10, linewidth=width, color=color, alpha=alpha)

    def show_plot(self):
        plt.show()

    def save_trajectory(self, file_name):
        plt.savefig(file_name, dpi=300)

    def start_animation(self):
        fig = self.fig
        ax = self.axes
        self.line, = ax.plot([], [], 'ro')
        self.ani = animation.FuncAnimation(fig, self._set_anim_pos, self._anim_data, blit=False, interval=10,
            repeat=True)
        plt.show()

    def save_animation(self, file_name):
        self.ani.save(file_name)

    def _anim_data(self):
        for xy in self.trajectory:
            yield xy[0]/10, xy[1]/10

    def _set_anim_pos(self, data):
        x, y = data[0], data[1]
        self.line.set_data(x, y)
        return self.line





