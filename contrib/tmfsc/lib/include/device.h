/** 
 * File: device.h
 * Author: K M Masum Habib
 * All dimensions assumed to be in nano-meters (nm);
 */

#ifndef TMFSC_LIB_DEVICE_H
#define TMFSC_LIB_DEVICE_H

#include "utils/Printable.hpp"
#include "maths/arma.hpp"
#include "segment.h"

#include <vector>

namespace qmicad{ namespace tmfsc {
static constexpr double nm = 1E-9;   // nanometer
static constexpr double AA = 1E-10;  //  Angstrom

using utils::Printable;
using maths::armadillo::linspace;
using maths::armadillo::vec;
using std::vector;

class Edge : public Segment {
public:
    Edge(const point &p, const point &q, int type = EDGE_REFLECT);

    void type(int type) { mType = type; }
    int type() { return mType; }

public:
    // edge types
    static constexpr int EDGE_REFLECT = 1;
    static constexpr int EDGE_ABSORB = 2;
    static constexpr int EDGE_TRANSMIT = 3;

protected:
    int mType;
};


class Device : public Printable {
public:
    Device();
    /** Adds point to the device. The device has a closed polygon, so add 
     *  first point again. */
    void addPoint(const point &pt);
    void addPoints(const vector<point> &pts);
    void edgeType(int iEdge, int type);
    int intersects(point p, point q);
    point intersection(int iEdge, point p, point q);
    vector<point> createPointsOnCont(int iCnt, int n);
    vector<point> createPointsOnCont(int iCnt, double dl);
    
    int contToEdgeIndx(int iCnt) { return mCnts[iCnt]; };
    int numEdges() { return mEdgs.size(); };
    int numConts() { return mCnts.size(); };
    double contDirctn(int iCnt) { return mEdgs[mCnts[iCnt]].angle(); }

    int edgeType(int iEdge) { return mEdgs[iEdge].type(); };
    bool isReflectEdge(int iEdge);
    bool isAbsorbEdge(int iEdge);
    const svec& edgeUnitVect(int indx) { return mEdgs[indx].unitVect(); };
    const svec edgeVect(int indx) { return mEdgs[indx].vect(); };
    const Edge& edge(int indx) { return mEdgs[indx]; };

public:
private:
    vector<point> mPts; //!< vertices.
    vector<Edge> mEdgs; //!< edges.
    vector<int> mCnts;  //!< contacts: index of absorbing edges.
};


}}

#endif


