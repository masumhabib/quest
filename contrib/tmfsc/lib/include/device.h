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
#include <map>

namespace qmicad{ namespace tmfsc {
using utils::Printable;
using maths::armadillo::linspace;
using maths::armadillo::vec;
using std::vector;
using std::map;

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
    int addPoint(const point &pt);
    void addPoints(const vector<point> &pts);
    void edgeType(int iEdge, int type);
    int intersects(point p, point q);
    point intersection(int iEdge, point p, point q);
    vector<point> createPointsOnCont(int iCnt, int n);
    vector<point> createPointsOnCont(int iCnt, double dl);
    
    int contToEdgeIndx(int iCnt) const { return mCnts[iCnt]; };
    int edgeToContIndx(int iEdge) const { return mEdg2Cnt.at(iEdge); };
    int numEdges() { return mEdgs.size(); };
    int numConts() { return mCnts.size(); };


    int edgeType(int iEdge) { return mEdgs[iEdge].type(); };
    bool isReflectEdge(int iEdge);
    bool isAbsorbEdge(int iEdge);
    const svec& edgeUnitVect(int indx) const { return mEdgs[indx].unitVect(); };
    svec edgeNormVect(int indx) const { return mEdgs[indx].normal(); };
    svec edgeVect(int indx) const { return mEdgs[indx].vect(); };
    const Edge& edge(int indx) const { return mEdgs[indx]; };
    
    
    svec contNormVect(int indx) const { 
        return edgeNormVect(contToEdgeIndx(indx)); };
    const svec& contUnitVect(int indx) const { 
        return edgeUnitVect(contToEdgeIndx(indx)); };
    double contDirctn(int iCnt) { return mEdgs[mCnts[iCnt]].angle(); }

public:
private:
    vector<point> mPts; //!< vertices.
    vector<Edge> mEdgs; //!< edges.
    vector<int> mCnts;  //!< Contact to edge map.
    map<int, int> mEdg2Cnt; //!< Egde to contact map.
};


}}

#endif


