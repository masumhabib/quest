/** file: segment.h
 * @author K M Masum Habib
 * Intersection algorithm of two line-segments.
 * Inspired by: 
 * http://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
"
 */

#ifndef TMFSC_LIB_SEGMENT_H
#define TMFSC_LIB_SEGMENT_H

#include "utils/Printable.hpp"
#include "maths/arma.hpp"
#include "maths/svec.h"
#include <stdexcept>


namespace qmicad{ namespace tmfsc{
using utils::Printable;
using maths::spvec::svec;
typedef maths::spvec::svec point;

using std::invalid_argument;
using std::max;
using std::min;
using std::abs;

class Segment : public Printable {
public:
    Segment(const point &p, const point &q);

    const svec& unitVect() const { return mn; };
    svec normal() const { return svec ({-mn(1), mn(0)}); };
    const svec vect() const { return svec(mr*mn); };
    double length() const { return mr; };
    double angle() const { return mth; };
    const point& p() const { return mp; };
    const point& q() const { return mq; };

    bool intersects(const point &p, const point &q, bool collinear = false);
    bool intersects(const Segment &seg, bool collinear = false);
    point intersection(const Segment &seg);
    point intersection(const point &p, const point &q);
 
protected:
    int orientation(const point &p, const point &q, const point &r);
    bool onSegment(const point &p, const point &q, const point &r);

protected:
    point mp;    // point 1
    point mq;    // point 2
    svec mn;     // unit vector pointing from q to p
    double mr;   // length of the vector
    double mth;  // direction of the vector
    double ma, mb, mc; // alternate representation

    static constexpr double TOL = 1E-10;
};

}}

#endif
