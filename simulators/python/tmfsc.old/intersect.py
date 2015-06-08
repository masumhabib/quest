"""
    FileName: intersect.py
    Intersection algorithm of two line-segments.
    Inspired by: http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
"""

import numpy as np

tol = 1E-10


def intersection(p1, q1, p2, q2):
    """
    Finds intersection point of two lines.
    Blatantly copied from SO: http://stackoverflow.com/questions/20677795/find-the-point-of-intersecting-lines
    """
    def line(p, q):
        a = (p[1] - q[1])
        b = (q[0] - p[0])
        c = p[0]*q[1] - q[0]*p[1]
        return a, b, -c

    L1 = line(p1, q1)
    L2 = line(p2, q2)
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if abs(D) > tol:
        x = Dx / D
        y = Dy / D
        return x, y
    else:
        return False


def intersects(p1, q1, p2, q2, collinear=True):
    """
        Returns true if two line segments, p1q1 and p2q2 intersect.
    """
    # orientations of all possible combinations of three points
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    # non-parallel and non-collinear case
    if o1 != o2 and o3 != o4:
        return True

    # For our case, we do not need to handle collinear
    if collinear:
        # collinear case
        if o1 == 0 and onsegment(p1, p2, q1):   # p2 lies on p1q1
            return True
        if o2 == 0 and onsegment(p1, q2, q1):   # q2 lies on p1q1
            return True
        if o3 == 0 and onsegment(p2, p1, q2):   # p1 lies on p2q2
            return True
        if o4 == 0 and onsegment(p2, q1, q2):   # q1 lies on p2q2
            return True

    return False


def orientation(p, q, r):
    """
        Finds orientation of ordered triplet(p,q,r)
        Returns:
            0 -- p,q,r are collinear
            1 -- Clockwise
            2 -- Counterclockwise
    """
    val = (q[1] - p[1])*(r[0]-q[0]) - (q[0]-p[0])*(r[1]-q[1])

    if abs(val) < tol:  # collinear case
        return 0
    if val > 0:         # clockwise
        return 1
    else:               # counterclockwise
        return 2


def onsegment(p, q, r):
    """
        Returns true if point q lies on segment pq, given three collinear
        points p, q, r.
    """
    if q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and q[1] <= max(p[1], r[1]) and q[1] >= min(p[1], r[1]):
        return True
    else:
        return False


def test(l1, l2):
    import matplotlib.pyplot as plt
    print ("Line 1:" + str(l1))
    print ("Line 2:" + str(l2))
    ints = intersects(l1[0,:], l1[1,:], l2[0,:], l2[1,:],False)
    if ints:
        print ("Intersection at: " + str(intersection(l1[0,:],l1[1,:],l2[0,:],l2[1,:])))
    else:
        print ("No intersection.")

    plt.plot(l2[:,0],l2[:,1])
    plt.show()
