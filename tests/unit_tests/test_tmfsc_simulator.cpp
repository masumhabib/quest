/** Test cases for Segment class.
 *
 */

#include "tmfsc/simulator.h"

#ifndef LINK_STATIC
#define BOOST_TEST_DYN_LINK
#endif
#define BOOST_TEST_MODULE SimulatorTest
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <memory>

using namespace quest::tmfsc;
using namespace std;

BOOST_AUTO_TEST_CASE(basics)
{

    point A = {0, 0};    
    point B = {20, 0};    
    point C = {20, 20};    
    point D = {0, 20};    

    Device::ptr dev = make_shared<Device>();
    dev->addPoint(A);
    dev->addPoint(B);
    dev->addPoint(C);
    dev->addPoint(D);

    dev->addEdge(0, 1);
    dev->addEdge(1, 2);
    dev->addEdge(2, 3);
    dev->addEdge(3, 0);

    dev->edgeType(1, Edge::EDGE_ABSORB);

    Simulator sim(dev);
    sim.setMaxNumStepsPerTraj(100);

    point p = {10, 5};
    TrajectoryVect trajs = sim.calcTraj(p, 0.0, 10, 0, 0.15);
    cout << " *** TRAJECTORY *** " << endl;
    for (auto traj : trajs) {
        for (int ip = 0; ip < traj.path.size(); ip += 1) {
            cout << "[" << ip << "]" << traj.path[ip] << endl;
        }
    }





    
}


