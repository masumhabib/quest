
#include "tmfsc/segment.h"
namespace quest { namespace tmfsc {

/** Creates a line segment joining (p,q) also provides an alternet 
 * representation in terms of co-efficients a,b and c */
Segment::Segment(const point &p, const point &q) 
    : mp(p), mq(q){
    ma = p[1] - q[1];
    mb = q[0] - p[0];
    mc = -(p[0]*q[1] - q[0]*p[1]);

    double dx = mb;
    double dy = -ma;
    mr = sqrt(dx*dx+dy*dy);
    mth = atan2(dy, dx);
    mn = {dx/mr, dy/mr};
}


/**
 * Returns true if two line segments, p1q1 and p2q2 intersect.
 */
bool Segment::intersects(const point &p, const point &q, bool collinear) {
   // orientations of all possible combinations of three points
    int o1 = orientation(mp, mq, p);
    int o2 = orientation(mp, mq, q);
    int o3 = orientation(p, q, mp);
    int o4 = orientation(p, q, mq);

    // non-parallel and non-collinear case
    if (o1 != o2 && o3 != o4) {
        return true;
    }

    // For our case, we do not need to handle collinear
    if (collinear) {
        // collinear case
        if (o1 == 0 && onSegment(mp, p, mq))   // p2 lies on p1q1
            return true;
        if (o2 == 0 && onSegment(mp, q, mq))   // q2 lies on p1q1
            return true;
        if (o3 == 0 && onSegment(p, mp, q))   // p1 lies on p2q2
            return true;
        if (o4 == 0 && onSegment(p, mq, q))   // q1 lies on p2q2
            return true;
    }

    return false;
}


bool Segment::intersects(const Segment &seg, bool collinear) {
    return intersects(seg.mp, seg.mq, collinear);
}

/**
 * Returns the intersection of two line segments. If intersection does 
 * not exist, throws error.
 * Blatantly copied from SO: 
 * http://stackoverflow.com/questions/20677795/find-the-point-of-intersecting-lines
 */
point Segment::intersection(const Segment &seg){
    double D  = ma * seg.mb - mb * seg.ma;
    double Dx = mc * seg.mb - mb * seg.mc;
    double Dy = ma * seg.mc - mc * seg.ma;
    if (abs(D) < TOL){
        throw invalid_argument("Lines L1 and L2 does not intersect");
    }

    point p = {Dx/D, Dy/D};
    return p;
}

point Segment::intersection(const point &p, const point &q) {
    return intersection(Segment(p, q));
}

/**
  * Finds orientation of ordered triplet(p,q,r)
  * Returns:
  *     0 -- p,q,r are collinear
  *     1 -- Clockwise
  *     2 -- Counterclockwise
  */
int Segment::orientation(const point &p, const point &q, const point &r) {
    double val = (q[1] - p[1])*(r[0]-q[0]) - (q[0]-p[0])*(r[1]-q[1]);
  
    // collinear case
    if (abs(val) < TOL) {
        return 0;
    }
    
    //clockwise
    if (val > 0) {
        return 1;
    }
    
    // counterclockwise
    return 2;
}

/**
  * Returns true if point q lies on segment pq, given three collinear
  * points p, q, r.
  */
bool Segment::onSegment(const point &p, const point &q, const point &r) {
    if (q[0] <= max(p[0], r[0])
           && q[0] >= min(p[0], r[0])
           && q[1] <= max(p[1], r[1]) 
           && q[1] >= min(p[1], r[1])) {
        return true;
    }

    return false;
}


}}
