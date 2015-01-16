
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
        self.pts = []
        self.edges = []
        self.edge_vec = []
        self.edge_types = []
        self.collected = []
        self.trajectory = np.array([])

    # dimensions are in angstroms
    def add_pts(self, pts):
        n = len(self.pts)
        i = n
        for pt in pts:
            self.pts.append(pt)
            self.collected.append(0)

            if i > 0:
                pt_im1 = self.pts[i-1]
                self.edges.append(np.array([pt_im1, pt]))
                dx = pt[0] - pt_im1[0]
                dy = pt[1] - pt_im1[1]
                r = sqrt(dx*dx+dy*dy)
                self.edge_vec.append(np.array([dx, dy])/r)
                self.edge_types.append(EDGE_REFLECT)

            i += 1

    def set_edge_type(self, indx, edge_type):
        self.edge_types[indx] = edge_type

    def intersects(self, pti, ptf):
        pti = pti/AA
        ptf = ptf/AA
        i = 0
        for edge in self.edges:
            if intersects(edge[0, :], edge[1,:], pti, ptf, False):
                pt = intersection(edge[0, :], edge[1,:], pti, ptf)
                if pt:
                    return i, pt[0]*AA, pt[1]*AA
            i += 1
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

    def draw_geom(self):
        self.fig = plt.figure()
        self.axes = self.fig.add_subplot(111)
        for ie in range(0, len(self.edges)):
            if self.edge_types[ie] == EDGE_ABSORB:
                width = 4.0
            else:
                width = 1.5
            edge = self.edges[ie]
            self.axes.plot(edge[:, 0]/10, edge[:, 1]/10, 'r-', linewidth=width)
        self.axes.set_aspect('equal', 'datalim')

    def draw_trajectory(self):
        width = 1.5
        self.axes.plot(self.trajectory[:, 0]/10, self.trajectory[:, 1]/10, linewidth=width)

    def show_plot(self):
        plt.show()

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
            # print xy[0]/10, xy[1]/10
            yield xy[0]/10, xy[1]/10

    def _set_anim_pos(self, data):
        x, y = data[0], data[1]
        self.line.set_data(x, y)
        return self.line





