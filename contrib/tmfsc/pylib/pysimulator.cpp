/** 
 * tmfsc::simulator python wrapper 
 */

#include "pysimulator.h"

namespace qmicad{ namespace python{

mat PySimulator::calcTrajPy(point ri, double thi, double B, 
            double EF, double V, bool saveTraj){
    auto points = calcTraj(ri, thi, B, EF, V, saveTraj);
    mat mpts(points.size(), 2);

    for (int itr = 0; itr < points.size(); ++itr){
        mpts.row(itr) = points[itr];
    }

    return mpts;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTrajPy, calcTrajPy, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTran, calcTran, 3, 4)
void export_Simulator(){
    using namespace qmicad::tmfsc;
    
    class_<PySimulator, bases<Printable>, shared_ptr<PySimulator> >("Simulator", 
            init<PyDevice&>())
        .def("calcTraj", &PySimulator::calcTrajPy, PySimulator_calcTrajPy())
        .def("calcTrans", &PySimulator::calcTran, PySimulator_calcTran())
    ;
}

}}


