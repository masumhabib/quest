/**
 * File: device.cpp
 * Author: K M Masum Habib
 */

#include "device.h"

namespace qmicad{ namespace tmfsc{

constexpr int Edge::EDGE_REFLECT;
constexpr int Edge::EDGE_ABSORB;
constexpr int Edge::EDGE_TRANSMIT;

/** Constructor */
Device::Device() {
}

/** Add a vertex to this device 
 */
void Device::addPoint(const point &pt) {
    mPts.push_back(pt);
    if (mPts.size() > 1){
        mEdgs.push_back(Edge(mPts[mPts.size() - 2], pt));
    }
}

/** Add vertices to this device. 
  */
void Device::addPoints(const vector<point> &pts) {
    for (point pt : pts){
        addPoint(pt);
    }
}

/** Add contact (edge with absorbing boundary)
 */
void Device::edgeType(int iEdge, int type) {
    mEdgs[iEdge].type(type);
    if (type == Edge::EDGE_ABSORB){
        mCnts.push_back(iEdge);
    }
}

/** Returns the edge index that intersects the line segment defined by
 *  p and q. p and q are in angstrom. */
int Device::intersects(point p, point q) {
    // convert to angstrom
    for (int i = 0; i < mEdgs.size(); i += 1) {
        if (mEdgs[i].intersects(p, q)) {
            return i;
        }
    }
    
    return -1;
}

/** Returns intersection point between and edge and the line segment defined 
 * by p and q. p and q are in angstrom. */
point Device::intersection(int iEdge, point p, point q) {
    // convert to angstrom
    return mEdgs[iEdge].intersection(p, q);
}

/** Returns N points on the contact contTndx.
 */
vector<point> Device::createPointsOnCont(int iCnt, int n) {
    Edge &edge = mEdgs[mCnts[iCnt]];
    point p = edge.p();
    point q = edge.q();
    vec x = linspace<vec>(p[0], q[0], n+2);
    vec y = linspace<vec>(p[1], q[1], n+2);

    vector<point> pts;
    for (int i = 1; i < n+1; i += 1){
        // shift the point a bit along the normal of the edge to avoid
        // registering on the injecting contact.
        // FIXME: probably this is not a good place for doing this.
        pts.push_back(svec({x[i], y[i]}) + edge.normal()/1E6);
    }

    return pts;
}

/** Returns N points on the contact contTndx.
 */
vector<point> Device::createPointsOnCont(int iCnt, double dl) {
    int npts = int(mEdgs[iCnt].length()/dl);
    return createPointsOnCont(iCnt, npts);
}
 
bool Device::isReflectEdge(int iEdge){
    return mEdgs[iEdge].type() == Edge::EDGE_REFLECT;
}

bool Device::isAbsorbEdge(int iEdge){
    return mEdgs[iEdge].type() == Edge::EDGE_ABSORB;
}


Edge::Edge(const point &p, const point &q, int type)
: Segment(p,q), mType(type){

}

}}




