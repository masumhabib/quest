/* 
 * File:   pyqmicad.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 3, 2014, 10:58 PM
 * 
 * Python interface to our wrappers.
 * 
 */


#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/wrap_python.hpp>

#include "../include/qmicad.hpp"

#include "pyqmicad.h"


namespace qmicad{
namespace python{

char const* greet()
{   
    static string msg;
    msg  = " QMICAD: Quantum Mechanics Inspired Computer Aided Design \n ";
    msg += "                        v" + qmicad::version;
    return msg.c_str();
}

//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Device_VDS, VDS, 1, 2)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Device_VLR, VLR, 3, 5)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyNegfEloop_enableTE, enableTE, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyGrapheneKpHam_setSize, setSize, 1, 2)
BOOST_PYTHON_MODULE(qmicad)
{
    using namespace boost::python;
    using boost::mpi::communicator;
    using boost::shared_ptr;
    using namespace std;
    
    def("greet", greet);
   
    
    /**
     * Our favorite armadillo matrices.
     */
    class_<cxmat, boost::shared_ptr<cxmat> >("cxmatp", boost::python::no_init)
    ;
    
    class_<mat, boost::shared_ptr<mat> >("matp", boost::python::no_init)
    ;

    
    /**
     * Periodic table.
     */
    class_<PyPeriodicTable, shared_ptr<PyPeriodicTable> >("PeriodicTable")
        .def("add", &PyPeriodicTable::add)
    ;

    /**
     * Atomistic geometry of the device.
     */
    class_<PyAtomicStruct, shared_ptr<PyAtomicStruct> >("AtomicStruct", init<const string&>())
        .def(init<const string&, const PyPeriodicTable>())
        .def(init<uint, uint, double, double, const PyPeriodicTable>())
        .def("span", &PyAtomicStruct::span)
        .def(self_ns::str(self_ns::self))
    ;
    
    /**
     * Graphene k.p parameters.
     */
    class_<PyGrapheneKpParams, shared_ptr<PyGrapheneKpParams> >("GrapheneKpParams")
        .def_readwrite("dtol", &PyGrapheneKpParams::dtol)
        .def_readwrite("ax", &PyGrapheneKpParams::ax)
        .def_readwrite("ay", &PyGrapheneKpParams::ay)
        .def_readwrite("K", &PyGrapheneKpParams::K)   
        .def_readwrite("gamma", &PyGrapheneKpParams::gamma)
        .def("update", &PyGrapheneKpParams::update)
        .def(self_ns::str(self_ns::self))
    ;

    /**
     * Graphene Hamiltonian.
     */
    //void (PyGrapheneKpHam::*PyGrapheneKpHam_generate1)(const AtomicStruct&, const AtomicStruct&, uint, uint) = &PyGrapheneKpHam::generate;
    class_<PyGrapheneKpHam, shared_ptr<PyGrapheneKpHam> >("GrapheneKpHam",init<const PyGrapheneKpParams& >())
        .def("setSize", &PyGrapheneKpHam::setSize, PyGrapheneKpHam_setSize())
        .def("getH0", &PyGrapheneKpHam::getH0)
        .def("getHl", &PyGrapheneKpHam::getHl)
        .def("getH", &PyGrapheneKpHam::getH)  
        .def("getSl", &PyGrapheneKpHam::getSl)
        .def("getS0", &PyGrapheneKpHam::getS0)
        .def("getS", &PyGrapheneKpHam::getS)  
        .def("genDiagBlock", &PyGrapheneKpHam::genDiagBlock)
        .def("genLowDiagBlock", &PyGrapheneKpHam::genLowDiagBlock)
        .def("generate", &PyGrapheneKpHam::generate)
    ;
    
    /*
    class_<DeviceParams>("DeviceParams")
        .def_readwrite("nl", &DeviceParams::nl)
        .def_readwrite("nw", &DeviceParams::nw)
        .def_readwrite("dtol", &DeviceParams::dtol)
        .def_readwrite("ax", &DeviceParams::ax)
        .def_readwrite("ay", &DeviceParams::ay)
        .def_readwrite("K", &DeviceParams::K)   
        .def_readwrite("gamma", &DeviceParams::gamma)
        .def_readwrite("A2", &DeviceParams::A2)
        .def_readwrite("C", &DeviceParams::C)
        .def_readwrite("dVD", &DeviceParams::dVDD)
        .def_readwrite("VDmax", &DeviceParams::VDDmx)
        .def_readwrite("VDmin", &DeviceParams::VDDmn)
        .def_readwrite("dVG", &DeviceParams::dVG)
        .def_readwrite("VGmax", &DeviceParams::VGmax)
        .def_readwrite("VGmin", &DeviceParams::VGmin)
        .def_readwrite("PotSolveType", &DeviceParams::PotSolveType)
        .def_readwrite("kT", &DeviceParams::kT)
        .def_readwrite("eta", &DeviceParams::eta)
        .def_readwrite("mu", &DeviceParams::mu)
        .def_readwrite("Emin", &DeviceParams::Emin)
        .def_readwrite("Emax", &DeviceParams::Emax)
        .def_readwrite("dE", &DeviceParams::dE)
        .def_readwrite("AutoGenE", &DeviceParams::AutoGenE)
        .def_readwrite("OutFileName", &DeviceParams::OutFileName)
        .def_readwrite("gjfFileName", &DeviceParams::gjfFileName)
        .def_readwrite("HamiltonianType", &DeviceParams::HamiltonianType)
        .def_readwrite("PotentialSolver", &DeviceParams::PotentialSolver)
        .def_readwrite("calcTE", &DeviceParams::CalcTE)
    ;

    class_<point>("Point", init<const double&, const double&>())
    ;
    
    class_<SimpleQuadrilateral>("Quadrilateral", init<const point&, const point&, const point&, const point&>())
    ;

    class_<Device>("Device", init<const communicator&, const DeviceParams&>())
        .def("tic", &Device::tic)
        .def("toc", &Device::toc)
        .def("time", &Device::time)
        .def("setDeviceParams", &Device::setDeviceParams)
        .def("prepare", &Device::prepare)
        .def("addSource", &Device::addSource)
        .def("addDrain", &Device::addDrain)
        .def("addGate", &Device::addGate)
        .def("addLinearRegion", &Device::addLinearRegion)
        .def("computePotential", &Device::computePotential)
        .def("runNegfEloop", &Device::runNegfEloop)
        .def("NegfParam", &Device::NegfParam)
        .def("toString", &Device::toString)
        .def("xmax", &Device::xmax)
        .def("xmin", &Device::xmin)
        .def("ymax", &Device::ymax)
        .def("ymin", &Device::ymin)
        .def("VDD", &Device::VDD)
        .def("NVDD", &Device::NVDD)
        .def("VGG", &Device::VGG)
        .def("NVGG", &Device::NVGG)
        .def("VDS", &Device::VDS, Device_VDS())
        .def("VG", &Device::VG)
        .def("VLR", &Device::VLR, Device_VLR())
    ; 
 
    class_<PyNegfEloop>("QnegfEloop", init<VecGrid&, const NegfParams&, const communicator&>())
        .def_readonly("TE", &PyNegfEloop::TE)
        .def_readonly("I1", &PyNegfEloop::I1)
        .def("prepare", &PyNegfEloop::prepare)
        .def("preCompute", &PyNegfEloop::preCompute)
        .def("compute", &PyNegfEloop::compute)
        .def("postCompute", &PyNegfEloop::postCompute)
        .def("collect", &PyNegfEloop::collect)
        .def("computeTE", &PyNegfEloop::computeTE)
        .def("collectTE", &PyNegfEloop::collectTE)
        .def("stepCompleted", &PyNegfEloop::stepCompleted)
        .def("enableTE", &PyNegfEloop::enableTE, PyNegfEloop_enableTE())
        .def("saveTE", &PyNegfEloop::saveTE)        
    ;

    class_<PyNegfParams>("QnegfParams", init<>())
    ;*/
    
}

}
}