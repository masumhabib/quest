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

tuple PySimulator::calcTranPy(double B, double E, double V, 
        int injCont, bool saveTraj){
    TrajectoryVect trajs;
    mat TE;
    tie(TE, trajs) = calcTran(B, E, V, injCont, saveTraj);

    list trajList; 
    int ntraj = trajs.size();
    for(int itraj = 0; itraj < ntraj; itraj += 1){
        int npts = trajs[itraj].size();
        mat pts (npts, 2);
        for (int ip = 0; ip < npts; ip += 1){
            pts.row(ip) = trajs[itraj][ip];
        }
        trajList.append(pts);
    }
    
    return make_tuple(TE, trajList);
}
 
int PySimulator::getParticleTypePy() {
    return static_cast<int>(getParticleType());
}

void PySimulator::setParticleTypePy(int type) {
    setParticleType(static_cast<ParticleType>(type));
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTrajPy, calcTrajPy, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTranPy, calcTranPy, 3, 5)
void export_Simulator(){
    using namespace qmicad::tmfsc;

    class_<PySimulator, bases<Printable>, shared_ptr<PySimulator> >("Simulator", 
            init<PyDevice&>())
        .def("calcTraj", &PySimulator::calcTrajPy, PySimulator_calcTrajPy())
        .def("calcTrans", &PySimulator::calcTranPy, PySimulator_calcTranPy())
        .add_property("ParticleType", &PySimulator::getParticleTypePy, 
                &PySimulator::setParticleTypePy)
        .add_property("NumPointsPerCycle", &PySimulator::getNumPointsPerCycle, 
                &PySimulator::getNumPointsPerCycle)

    ;
}

}}


