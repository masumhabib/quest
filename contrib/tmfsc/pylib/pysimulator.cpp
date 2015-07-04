/** 
 * tmfsc::simulator python wrapper 
 */

#include "pysimulator.h"

namespace qmicad{ namespace python{

mat PySimulator::calcTrajPy(point ri, double thi, double B, 
            double EF, double V, bool saveTraj)
{
    auto points = calcTraj(ri, thi, B, EF, V, saveTraj);
    return Traj2Mat(points);
}

tuple PySimulator::calcTranPy(double B, double E, double V, 
        int injCont, bool saveTraj)
{
    TrajectoryVect trajs;
    mat TE;
    tie(TE, trajs) = calcTran(B, E, V, injCont, saveTraj);

    list trajList; 
    int ntraj = trajs.size();
    for(int itraj = 0; itraj < ntraj; itraj += 1){
        trajList.append(Traj2Mat(trajs[itraj]));
    }
    
    return make_tuple(TE, trajList);
}
 
mat PySimulator::calcTrajPy2(point ri, double thi, double B, double E, 
        const list& VG, bool saveTraj) 
{
    vector<double> VGs = list2vect<double> (VG);
    auto points = calcTraj(ri, thi, B, E, VGs, saveTraj);

    return Traj2Mat(points);

}

tuple PySimulator::calcTranPy2(double B, double E, const list& VG, 
        int injCont, bool saveTraj) 
{
    TrajectoryVect trajs;
    mat TE;
    vector<double> VGs = list2vect<double> (VG);
    tie(TE, trajs) = calcTran(B, E, VGs, injCont, saveTraj);

    list trajList; 
    int ntraj = trajs.size();
    for(int itraj = 0; itraj < ntraj; itraj += 1){
        trajList.append(Traj2Mat(trajs[itraj]));
    }
    
    return make_tuple(TE, trajList);

}

int PySimulator::getParticleTypePy() {
    return static_cast<int>(getParticleType());
}

void PySimulator::setParticleTypePy(int type) {
    setParticleType(static_cast<ParticleType>(type));
}

mat PySimulator::Traj2Mat (const Trajectory& traj) {
    mat m(traj.size(), 2);
    for (int itr = 0; itr < traj.size(); ++itr){
        m.row(itr) = traj[itr];
    }
    return m;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTrajPy, calcTrajPy, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTranPy, calcTranPy, 3, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTrajPy2, calcTrajPy2, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTranPy2, calcTranPy2, 3, 5)
void export_Simulator(){
    using namespace qmicad::tmfsc;

    class_<PySimulator, bases<Printable>, shared_ptr<PySimulator> >("Simulator", 
            init<shared_ptr<PyDevice> >())
        .def("calcTraj", &PySimulator::calcTrajPy, PySimulator_calcTrajPy())
        .def("calcTrans", &PySimulator::calcTranPy, PySimulator_calcTranPy())
        .def("calcTraj", &PySimulator::calcTrajPy2, PySimulator_calcTrajPy2())
        .def("calcTrans", &PySimulator::calcTranPy2, PySimulator_calcTranPy2())
        .add_property("ParticleType", &PySimulator::getParticleTypePy, 
                &PySimulator::setParticleTypePy)
        .add_property("NumPointsPerCycle", &PySimulator::getNumPointsPerCycle, 
                &PySimulator::setNumPointsPerCycle)
        .add_property("TimeStep", &PySimulator::getTimeStep, 
                &PySimulator::setTimeStep)
    ;
}

}}


