/**
 * File: device.cpp
 * Author: K M Masum Habib
 * Co-Author: Mirza Elahi
 */

#include "device.h"

namespace qmicad{ namespace tmfsc{

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

/** Adds an edge devined by two point indices. */
int Device::addEdge(int ipt1, int ipt2, int type) {
    if (ipt1 >= mPts.size() || ipt2 >= mPts.size()) {
        throw invalid_argument(" Device::addEdge(): indices out of bounds");
    }
    mEdgs.push_back(Edge(mPts[ipt1], mPts[ipt2], type));
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

tuple<double, double, double, double> Device::calcProbab(double V1, double V2,
            const svec& vel, double En, int iEdge) 
{
	using maths::constants::e;
	using maths::constants::hbar;
	double th, thti, TransProb, RefProb;
	svec NormVec = this->edgeNormVect(iEdge);
	double incidentAbsoluteAngle = atan2( vel[1], vel[0] );
	//TODO have to think about a generalized formula
	if(incidentAbsoluteAngle<=pi && incidentAbsoluteAngle>=-pi/2){
		thti = incidentAbsoluteAngle - atan2( NormVec[1], NormVec[0]);
		if( abs(thti) > pi/2 ){
			thti += pi;
		}
	}
	// correcting theta incidence due to vector complication
	else{
		thti = atan2( NormVec[1], NormVec[0]) + incidentAbsoluteAngle;
	}

	// handle the critical angle case
	double sininv = abs((En+V1)/(En+V2)) * sin(thti);
	if (abs(sininv) > 1) {
	    return make_tuple(std::nan("NaN"), thti, 0.0, 1.0);
    }

	th = asin(abs((En+V1)/(En+V2)) * sin(thti));
    // if pn case, calculate exponential term and apply negative refraction
	double t_graded = 1.0;
	if ( !( (En > -V1 && En < -V2) || (En > -V2 && En < -V1) ) ){
		th = -th;
	} else {
	    //FIXME: REMOVE hardcoded vF, else, it is gonna cause major headache
    	double vF = 1E6;
    	double kf1, ky, d_eff, S;
    	kf1 = e * abs(En - (-V1)) / (hbar * vF);
    	ky = kf1 * sin ( thti );
    	d_eff = hbar * vF * abs(ky) * (this->splitLen*nm) / ( e*abs(V1-V2) );
    	S = pi * (d_eff/2) * abs(ky);
    	t_graded = std::exp( -2*S );
    }

	TransProb = t_graded*std::cos( thti )*std::cos( th ) 
	    /   std::pow(  cos( (abs(thti)+abs(th))/2 ), 2  )  ;

//    if( En+V2 > En+V1 )
//    {
//    	//TODO
//
//    }
	RefProb = 1 - TransProb;
	return make_tuple(th, thti, TransProb, RefProb);
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

Edge::Edge(const point &p, const point &q, int type)
: Segment(p,q), mType(type){

}

double Device::getSplitLen(){
	return this->splitLen;
}

void Device::setSplitLen(double len){
	this->splitLen = len;
}

}}




