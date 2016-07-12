/** Test cases for Segment class.
 *
 */

#include "segment.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SegmentTest
#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace qmicad::tmfsc;
using namespace std;
BOOST_AUTO_TEST_CASE(basic)
{
    point p = {0, 0};
    point q = {10, 10};
    Segment pq(p, q);

    BOOST_CHECK(pq.length() == sqrt(200));
    BOOST_CHECK_CLOSE(pq.angle(), 3.1415926536/4, 1E-4);

    svec unitvec = pq.unitVect();
    BOOST_CHECK_CLOSE(unitvec(0), 1/sqrt(2), 1E-3);
    BOOST_CHECK_CLOSE(unitvec(1), 1/sqrt(2), 1E-3);

    svec vect = pq.vect();
    BOOST_CHECK_CLOSE(vect(0), 10, 1E-3);
    BOOST_CHECK_CLOSE(vect(1), 10, 1E-3);
}

BOOST_AUTO_TEST_CASE(intersection)
{
    point p1 = {0, 0};
    point q1 = {10, 10};
    Segment seg1(p1, q1);

    point p2 = {0, 10};    
    point q2 = {10, 0};    

    BOOST_CHECK(seg1.intersects(p2, q2));

    point ints = seg1.intersection(p2, q2);
    BOOST_CHECK_CLOSE(ints(0), 5, 1E-4);
    BOOST_CHECK_CLOSE(ints(1), 5, 1E-4);

}



