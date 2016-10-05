/** 
 * tmfsc::simulator python wrapper 
 */

#include "pysimulator.h"

namespace quest{ namespace python{

list PySimulator::calcTrajPy(point ri, double thi, double E, 
            double B, double V, bool saveTraj)
{
    TrajectoryVect trajs = calcTraj(ri, thi, E, B, V, saveTraj);
    return TrajVect2List(trajs);
}

list PySimulator::calcTrajPy2(point ri, double thi, double E, double B, 
        const list& VG, bool saveTraj) 
{
    vector<double> VGs = list2vect<double> (VG);
    TrajectoryVect trajs = calcTraj(ri, thi, E, B, VGs, saveTraj);

    return TrajVect2List(trajs);

}

tuple PySimulator::calcTranPy(double E, double B, double V, 
        int injCont, bool saveTraj)
{
    TrajectoryVect trajs;
    mat TE;
    tie(TE, trajs) = calcTran(E, B, V, injCont, saveTraj);

    return make_tuple(TE, TrajVect2List(trajs));
}
 

tuple PySimulator::calcTranPy2(double E, double B, const list& VG, 
        int injCont, bool saveTraj) 
{
    TrajectoryVect trajs;
    mat TE;
    vector<double> VGs = list2vect<double> (VG);
    tie(TE, trajs) = calcTran(E, B, VGs, injCont, saveTraj);

    return make_tuple(TE, TrajVect2List(trajs));
}

int PySimulator::getParticleTypePy() {
    return static_cast<int>(getParticleType());
}

void PySimulator::setParticleTypePy(int type) {
    setParticleType(static_cast<ParticleType>(type));
}

int PySimulator::getInjectModelPy() {
    return static_cast<int>(getInjectModel());
}

void PySimulator::setInjectModelPy(int model) {
    setInjectModel(static_cast<InjectModel>(model));
}

list PySimulator::TrajVect2List(const TrajectoryVect& trajs) {
    list trajList; 
    int ntraj = trajs.size();
    for(int itraj = 0; itraj < ntraj; itraj += 1){
        trajList.append(Traj2PyTraj(trajs[itraj]));
    }

    return trajList;
}

PyTrajectory PySimulator::Traj2PyTraj (const Trajectory& traj) {
    PyTrajectory pytraj;
    pytraj.path.set_size(traj.path.size(),2);
    for (int itr = 0; itr < traj.path.size(); ++itr){
        pytraj.path.row(itr) = traj.path[itr];
    }
    pytraj.occupation = traj.occupation;
    return pytraj;
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTrajPy, calcTrajPy, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTranPy, calcTranPy, 3, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTrajPy2, calcTrajPy2, 5, 6)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PySimulator_calcTranPy2, calcTranPy2, 3, 5)
void export_Simulator(){
    using namespace quest::tmfsc;

    class_<PyTrajectory, shared_ptr<PyTrajectory> >("Trajectory")
        .add_property("path", &PyTrajectory::getPath)
        .def_readwrite("occupation", &PyTrajectory::occupation)
    ;

    class_<PySimulator, bases<Printable>, shared_ptr<PySimulator>, noncopyable>("Simulator", 
            init<shared_ptr<PyDevice> >())
        .def("calcTraj", &PySimulator::calcTrajPy, PySimulator_calcTrajPy())
        .def("calcTrans", &PySimulator::calcTranPy, PySimulator_calcTranPy())
        .def("calcTraj", &PySimulator::calcTrajPy2, PySimulator_calcTrajPy2())
        .def("calcTrans", &PySimulator::calcTranPy2, PySimulator_calcTranPy2())
        .add_property("MaxNumStepsPerTraj", &PySimulator::getMaxNumStepsPerTraj, 
                &PySimulator::setMaxNumStepsPerTraj)
        .add_property("NumPointsPerCycle", &PySimulator::getNumPointsPerCycle, 
                &PySimulator::setNumPointsPerCycle)
        .add_property("NumEdgeSeekSteps", &PySimulator::getNumEdgeSeekSteps, 
                &PySimulator::setNumEdgeSeekSteps)
        .add_property("InjecSpacing", &PySimulator::getInjecSpacing, 
                &PySimulator::setInjecSpacing)
        .add_property("NumInjecDir", &PySimulator::getNumInjecDir, 
                &PySimulator::setNumInjecDir)
        .add_property("InjecAngleSpread", &PySimulator::getInjecAngleSpread, 
                &PySimulator::setInjecAngleSpread)
        .add_property("MeanInjecAngle", &PySimulator::getMeanInjecAngle, 
                &PySimulator::setMeanInjecAngle)
        .add_property("ParticleType", &PySimulator::getParticleTypePy, 
                &PySimulator::setParticleTypePy)
        .add_property("TimeStep", &PySimulator::getTimeStep, 
                &PySimulator::setTimeStep)
        .add_property("CollectionTol", &PySimulator::getCollectionTol, 
                &PySimulator::setCollectionTol)
        .add_property("DebugLevel", &PySimulator::getDebugLvl,
                &PySimulator::setDebugLvl)
        .add_property("MinmNoInjection", &PySimulator::getMinNoInjection,
                  &PySimulator::setMinNoInjection)
        .add_property("EdgRghAngleSpread", &PySimulator::getEdgRghAngleSpread,
                  &PySimulator::setEdgRghAngleSpread)

    ;
}

}}


