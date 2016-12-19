/** Test cases for Segment class.
 *
 */

#include "tmfsc/device.h"

#ifndef LINK_STATIC
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE DeviceTest
#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace quest::tmfsc;
using namespace std;

Device create_simple_device (const point& A, const point& B, const point& C,
const point& D) {
    Device dev;
    int ptA = dev.addPoint(A);
    int ptB = dev.addPoint(B);
    int ptC = dev.addPoint(C);
    int ptD = dev.addPoint(D);

    dev.addEdge(ptA, ptB);
    dev.addEdge(ptB, ptC);
    dev.addEdge(ptC, ptD);
    dev.addEdge(ptD, ptA);

    return dev;
}

BOOST_AUTO_TEST_CASE(basics)
{
    Device dev = create_simple_device ({0, 0}, {20, 0}, {20, 10}, {0, 10});
    dev.edgeType(1, Edge::EDGE_ABSORB);
    
    BOOST_CHECK(dev.numEdges() == 4);
    BOOST_CHECK(dev.numConts() == 1);
    BOOST_CHECK(dev.isAbsorbEdge(1));
    BOOST_CHECK(dev.isReflectEdge(2));
    BOOST_CHECK(dev.contToEdgeIndx(0) == 1);
}

BOOST_AUTO_TEST_CASE(vectors)
{
    Device dev = create_simple_device ({0, 0}, {20, 0}, {20, 10}, {0, 10});
    dev.edgeType(1, Edge::EDGE_ABSORB);
    
    svec u = dev.edgeUnitVect(0);
    BOOST_CHECK_CLOSE(u(0), 1, 1E-4);
    BOOST_CHECK_CLOSE(u(1), 0, 1E-4);

    u = dev.edgeUnitVect(1);
    BOOST_CHECK_CLOSE(u(0), 0, 1E-4);
    BOOST_CHECK_CLOSE(u(1), 1, 1E-4);

    u = dev.edgeUnitVect(2);
    BOOST_CHECK_CLOSE(u(0), -1, 1E-4);
    BOOST_CHECK_CLOSE(u(1), 0, 1E-4);

    u = dev.edgeUnitVect(3);
    BOOST_CHECK_CLOSE(u(0), 0, 1E-4);
    BOOST_CHECK_CLOSE(u(1), -1, 1E-4);

    u = dev.edgeVect(3);
    BOOST_CHECK_CLOSE(u(0), 0, 1E-4);
    BOOST_CHECK_CLOSE(u(1), -10, 1E-4);
}

BOOST_AUTO_TEST_CASE(intersection)
{
    Device dev = create_simple_device ({0, 0}, {20, 0}, {20, 10}, {0, 10});
    dev.edgeType(1, Edge::EDGE_ABSORB);
    
    point ctr = {10, 5};
    point q = {30, 5};

    BOOST_CHECK(dev.intersects(ctr, q) == 1);
    point ints = dev.intersection(1, ctr, q);
    BOOST_CHECK_CLOSE(ints(0), 20, 1E-4);
    BOOST_CHECK_CLOSE(ints(1), 5, 1E-4);

    q(0) = 20; q(1) = 5;
    BOOST_CHECK(dev.intersects(ctr, q) == 1);

    q(0) = 20-0.0001; q(1) = 5;
    BOOST_CHECK(dev.intersects(ctr, q) == -1);

    point p1 = {20-0.01, 10-0.01};
    point q1 = {20-0.01, 10+0.01};
    BOOST_CHECK(dev.intersects(p1, q1) == 2);
    ints = dev.intersection(2, p1, q1);
    BOOST_CHECK_CLOSE(ints(0), 20-0.01, 1E-4);
    BOOST_CHECK_CLOSE(ints(1), 10, 1E-4);
}

BOOST_AUTO_TEST_CASE(pointsOnEdge)
{
    Device dev = create_simple_device ({0, 0}, {20, 0}, {20, 10}, {0, 10});
    dev.edgeType(1, Edge::EDGE_ABSORB);

    Edge e = dev.edge(1);
    cout << "Edge#0 p=(" << e.p()[0] << "," << e.p()[1] << ")" 
                 << " q=(" << e.q()[0] << "," << e.q()[1] << ")" 
                 << endl << endl; 

    vector<point> pts = dev.createPointsOnCont(0, 11);
    for (int i = 0; i < pts.size(); i += 1) {
        cout << "pt#" << i << pts[i];
    }
}




