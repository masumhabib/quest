/** Test cases for Segment class.
 *
 */

#include "device.h"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE DeviceTest
#include <boost/test/unit_test.hpp>

#include <iostream>

using namespace qmicad::tmfsc;
using namespace std;

BOOST_AUTO_TEST_CASE(basics)
{

    point A = {0, 0};    
    point B = {20, 0};    
    point C = {20, 10};    
    point D = {0, 10};    

    Device dev;
    dev.addPoint(A);
    dev.addPoint(B);
    dev.addPoint(C);
    dev.addPoint(D);
    dev.addPoint(A);

    dev.edgeType(1, Edge::EDGE_ABSORB);
    
    BOOST_CHECK(dev.numEdges() == 4);
    BOOST_CHECK(dev.numConts() == 1);
    BOOST_CHECK(dev.isAbsorbEdge(1));
    BOOST_CHECK(dev.isReflectEdge(2));
    BOOST_CHECK(dev.contToEdgeIndx(0) == 1);
}

BOOST_AUTO_TEST_CASE(vectors)
{

    point A = {0, 0};    
    point B = {20, 0};    
    point C = {20, 10};    
    point D = {0, 10};    

    Device dev;
    dev.addPoint(A);
    dev.addPoint(B);
    dev.addPoint(C);
    dev.addPoint(D);
    dev.addPoint(A);

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
    point A = {0, 0};    
    point B = {20, 0};    
    point C = {20, 10};    
    point D = {0, 10};    

    Device dev;
    dev.addPoint(A);
    dev.addPoint(B);
    dev.addPoint(C);
    dev.addPoint(D);
    dev.addPoint(A);

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

    point A = {0, 0};    
    point B = {20, 0};    
    point C = {20, 10};    
    point D = {0, 10};    

    Device dev;
    dev.addPoint(A);
    dev.addPoint(B);
    dev.addPoint(C);
    dev.addPoint(D);
    dev.addPoint(A);
    
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




