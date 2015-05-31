/** 
 * tmfsc::simulator python wrapper 
 */

#include "pysimulator.h"

namespace qmicad{ namespace python{

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Simulator_calcTraj, calcTraj, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Simulator_calcTran, calcTran, 3, 4)
void export_Simulator(){
    using namespace qmicad::tmfsc;
    class_<Simulator, bases<Printable>, shared_ptr<Simulator> >("Simulator", 
            init<Device&>())
        .def("calcTraj", &Simulator::calcTraj, Simulator_calcTraj())
        .def("calcTrans", &Simulator::calcTran, Simulator_calcTran())
    ;
}

}}


