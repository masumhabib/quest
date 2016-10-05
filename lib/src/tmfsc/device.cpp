/**
 * File: device.cpp
 * Author: K M Masum Habib
 * Co-Author: Mirza Elahi
 */

#include "tmfsc/device.h"

namespace quest{ namespace tmfsc{

constexpr int Edge::EDGE_REFLECT;
constexpr int Edge::EDGE_ABSORB;
constexpr int Edge::EDGE_TRANSMIT;

/** Constructor */
Device::Device() {
    mPot = make_shared<LinearPot>();
}

/** Add a vertex to this device 
 */
int Device::addPoint(const point &pt) {
    mPts.push_back(pt);
    /*if (mPts.size() > 1){
        mEdgs.push_back(Edge(mPts[mPts.size() - 2], pt));
        return mEdgs.size() - 1;
    }*/
    return mPts.size() - 1;
}

/** Add vertices to this device. 
  */
void Device::addPoints(const vector<point> &pts) {
    for (point pt : pts){
        addPoint(pt);
    }
}

/** Adds an edge defined by two point indices. */
int Device::addEdge(int ipt1, int ipt2, int type, double d) {
    if (ipt1 >= mPts.size() || ipt2 >= mPts.size()) {
        throw invalid_argument(" Device::addEdge(): indices out of bounds");
    }
    if (d < 0){
    	throw invalid_argument("Split length of junction cannot be negative");
    }
    if (type != Edge::EDGE_TRANSMIT && d>0){
    	throw invalid_argument("Edge other than transmitting cannot have split length");
    }
    mEdgs.push_back(Edge(mPts[ipt1], mPts[ipt2], type, d));
    return mEdgs.size() - 1;
}

/** Add contact (edge with absorbing boundary)
 */
void Device::edgeType(int iEdge, int type) {
    mEdgs[iEdge].type(type);
    if (type == Edge::EDGE_ABSORB){
        mCnts.push_back(iEdge);
        mEdg2Cnt[iEdge] = mCnts.size()-1;
    }
}

/** Returns the edge index that intersects the line segment defined by
 *  p and q. */
int Device::intersects(const point &p, const point &q) {
    for (int i = 0; i < mEdgs.size(); i += 1) {
        if (intersects(i, p, q)) {
            return i;
        }
    }
    
    return -1;
}

/** Returns true if edge iEdge intersects the segment pq. */
int Device::intersects(int iEdge, const point &p, const point &q) {
    return mEdgs[iEdge].intersects(p, q);
}

/** Returns intersection point between and edge and the line segment defined 
 * by p and q. p and q are in angstrom. */
point Device::intersection(int iEdge, const point &p, const point &q) {
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
    int npts = int(mEdgs[mCnts[iCnt]].length()/dl);
    return createPointsOnCont(iCnt, npts);
}
 
bool Device::isReflectEdge(int iEdge) {
    return mEdgs[iEdge].type() == Edge::EDGE_REFLECT;
}

bool Device::isAbsorbEdge(int iEdge) {
    return mEdgs[iEdge].type() == Edge::EDGE_ABSORB;
}

bool Device::isTransmitEdge(int iEdge) {
    return mEdgs[iEdge].type() == Edge::EDGE_TRANSMIT;
}

/** Transmission probability of an electron going through a transmitting
 *  edge  */
double Device::calcTransProb(double V1, double V2, const svec& vel, double En,
		int iEdge)
{
	using maths::constants::e;
	using maths::constants::hbar;
	double thtr, thti, TransProb;
	svec NormVec = this->edgeNormVect(iEdge);
	// incident angle
	double vMag = sqrt( vel[0]*vel[0] + vel[1]*vel[1] );
	svec v_unit = vel / vMag;
	thti = acos( dot(NormVec, v_unit) );
	// correction for incident angle for normal direction not in traditional
	// direction
	if( thti > pi/2 ){
		thti = pi - thti;
	}else{
		NormVec = -NormVec;
	}
	// total internal reflection
	double sininv = abs((En+V1)/(En+V2)) * sin(thti);
	if ( abs(sininv) > 1 ) {
		return 0.0;
	}
	// refraction angle
	thtr = asin(  ( (En+V1)/(En+V2) )  *  sin(thti)  );
	// if pn case, calculate exponential term
	double t_graded = 1.0;
	if (  (En > -V1 && En < -V2) || (En > -V2 && En < -V1) ){
		//FIXME: REMOVE hardcoded vF, else, it is gonna cause major headache
		double vF = 1E6;
		double kf1, ky, d_eff, S;
		kf1 = e * abs(En - (-V1)) / (hbar * vF);
		ky = kf1 * sin ( thti );
		d_eff = hbar * vF * abs(ky) * (mEdgs[iEdge].split()*nm) / ( e*abs(V1-V2) );
		S = pi * (d_eff/2) * abs(ky);
		t_graded = std::exp( -2*S );
	}
	// calc trans. prob.
	TransProb = t_graded  *  cos(thti) * cos(thtr)
				/ pow(  cos( (abs(thti)+abs(thtr))/2 ), 2  )  ;
	return TransProb;
}
 
int Device::addGate(const point& lb, const point& rb, const point& rt, 
        const point& lt) 
{
    typedef maths::geometry::point gpoint;
    return mPot->addGate(squadrilateral(gpoint(lb[0], lb[1]), 
                gpoint(rb[0], rb[1]), gpoint(rt[0], rt[1]), 
                gpoint(lt[0], lt[1])));
}

void Device::setGatePotential(int igate, double V) {
    mPot->VG(igate, V);
}

double Device::getPotAt(const point& position) {
    typedef maths::geometry::point gpoint;
    return mPot->getPotAt(gpoint(position[0], position[1]));
}

Edge::Edge(const point &p, const point &q, int type, double d)
: Segment(p,q), mType(type), md(d){

}
}}




