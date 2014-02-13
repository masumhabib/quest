/* 
 * File:   pyqmicad.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 3, 2014, 10:58 PM
 * 
 * Python interface to our wrappers.
 * 
 */

#include "pyqmicad.h"
#include "../include/qmicad.hpp"

#include <python2.6/Python.h>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/detail/wrap_python.hpp>



namespace qmicad{
namespace python{

char const* greet()
{   
    static string msg;
    msg  = "      QMICAD: Quantum Mechanics Inspired Computer Aided Design";
    msg += " v" + qmicad::version + "\n";
    msg += " -----------------------------------------------------------------------";
    return msg.c_str();
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyLinearPot_VLR, VLR, 3, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyNegfEloop_enableTE, enableTE, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyGrapheneKpHam_setSize, setSize, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyTISurfHam_setSize, setSize, 1, 2)
BOOST_PYTHON_MODULE(qmicad)
{
    using namespace boost::python;
    using boost::mpi::communicator;
    using boost::shared_ptr;
    using boost::noncopyable;
    using namespace utils::stds;
    
    def("greet", greet);
   
    /**
     * Enums.
     */
    enum_<Option>("Option") 
        .value("Disabled", Disabled)
        .value("Enabled",  Enabled)
    ;
    
    /**
     * Our favorite armadillo matrices.
     */
    class_<cxmat, shared_ptr<cxmat> >("cxmatp", no_init)
    ;
    
    class_<mat, shared_ptr<mat> >("matp", no_init)
    ;

    class_<vec, shared_ptr<vec> >("vecp", no_init)
    ;

    /**
     * Geometry classes
     */
    class_<point>("Point", init<const double&, const double&>())
    ;
    
    class_<SimpleQuadrilateral>("Quadrilateral", init<const point&, const point&, 
            const point&, const point&>())
    ;

    /**
     * Wall clock
     */
    class_<PyTimer, shared_ptr<PyTimer> >("Timer")
        .def("tic", &PyTimer::tic)
        .def("toc", &PyTimer::toc)
        .def(self_ns::str(self_ns::self))        
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
    class_<PyAtomicStruct, shared_ptr<PyAtomicStruct> >("AtomicStruct", 
            init<const string&>())
        .def(init<const string&, const PyPeriodicTable>())
        .def(init<uint, uint, double, double, const PyPeriodicTable>())
        .def("span", &PyAtomicStruct::span)
        .def("xmax", &PyAtomicStruct::xmax)
        .def("xmin", &PyAtomicStruct::xmin)
        .def("ymax", &PyAtomicStruct::ymax)
        .def("ymin", &PyAtomicStruct::ymin)    
        .def("zmax", &PyAtomicStruct::zmax)
        .def("zmin", &PyAtomicStruct::zmin)    
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
    class_<PyGrapheneKpHam, shared_ptr<PyGrapheneKpHam> >("GrapheneKpHam",
            init<const PyGrapheneKpParams& >())
        .def("setSize", &PyGrapheneKpHam::setSize, PyGrapheneKpHam_setSize())
        .def("H0", &PyGrapheneKpHam::H0)
        .def("Hl", &PyGrapheneKpHam::Hl)
        .def("H", &PyGrapheneKpHam::H)  
        .def("Sl", &PyGrapheneKpHam::Sl)
        .def("S0", &PyGrapheneKpHam::S0)
        .def("S", &PyGrapheneKpHam::S)  
        .def("genDiagBlock", &PyGrapheneKpHam::genDiagBlock)
        .def("genLowDiagBlock", &PyGrapheneKpHam::genLowDiagBlock)
        .def("generate", &PyGrapheneKpHam::generate)
    ;
    
    /**
     * TI Surface k.p parameters.
     */
    class_<PyTISurfKpParams, shared_ptr<PyTISurfKpParams> >("TISurfKpParams")
        .def_readwrite("dtol", &PyTISurfKpParams::dtol)
        .def_readwrite("ax", &PyTISurfKpParams::ax)
        .def_readwrite("ay", &PyTISurfKpParams::ay)
        .def_readwrite("K", &PyTISurfKpParams::K)   
        .def_readwrite("A2", &PyTISurfKpParams::A2)
        .def_readwrite("C", &PyTISurfKpParams::C)
        .def("update", &PyTISurfKpParams::update)
        .def(self_ns::str(self_ns::self))
    ;

    /**
     * TI Surface Hamiltonian.
     */
    class_<PyTISurfHam, shared_ptr<PyTISurfHam> >("TISurfHam",
            init<const PyTISurfHam& >())
        .def("setSize", &PyTISurfHam::setSize, PyTISurfHam_setSize())
        .def("H0", &PyTISurfHam::H0)
        .def("Hl", &PyTISurfHam::Hl)
        .def("H", &PyTISurfHam::H)  
        .def("Sl", &PyTISurfHam::Sl)
        .def("S0", &PyTISurfHam::S0)
        .def("S", &PyTISurfHam::S)  
        .def("genDiagBlock", &PyTISurfHam::genDiagBlock)
        .def("genLowDiagBlock", &PyTISurfHam::genLowDiagBlock)
        .def("generate", &PyTISurfHam::generate)
    ;
    
    /**
     * Linear potential
     */    
    class_<PyLinearPot, shared_ptr<PyLinearPot> >("LinearPot", 
            init<const PyAtomicStruct&, optional<const string&> >())
        .def("addSource", &PyLinearPot::addSource)
        .def("addDrain", &PyLinearPot::addDrain)
        .def("addGate", &PyLinearPot::addGate)
        .def("addLinearRegion", &PyLinearPot::addLinearRegion)
        .def("compute", &PyLinearPot::compute)
        .def("exportSvg", &PyLinearPot::exportSvg)
        .def("exportPotential", &PyLinearPot::exportPotential)
        .def("VD", &PyLinearPot::VD)
        .def("VS", &PyLinearPot::VS)
        .def("VG", &PyLinearPot::VG)
        .def("VLR", &PyLinearPot::VLR, PyLinearPot_VLR()) 
        .def("toOrbPot", &PyLinearPot::toOrbPot) 
        .def(self_ns::str(self_ns::self))
    ;
     
    /**
     * Vector grid
     */
    class_<PyVecGrid, shared_ptr<PyVecGrid> >("VecGrid", 
            init<double, double, double, optional<const string&> >())
        .def(init<optional<double, double, int, const string&> >())
        .def("V", &PyVecGrid::V)        
        .def("min", &PyVecGrid::min)
        .def("max", &PyVecGrid::max)
        .def("del", &PyVecGrid::del)
        .def("N", &PyVecGrid::N)
        .def(self_ns::str(self_ns::self))
    ;

    /**
     * NEGF parameters.
     */
    class_<PyNegfParams, shared_ptr<PyNegfParams> >("NegfParams", init<uint>())
        .def("H0", &PyNegfParams::H0)
        .def("S0", &PyNegfParams::S0)
        .def("Hl", &PyNegfParams::Hl)
        .def("Sl", &PyNegfParams::Sl)
        .def("V", &PyNegfParams::V)
        .def_readwrite("kT", &PyNegfParams::kT)
        .def_readwrite("ieta", &PyNegfParams::ieta)
        .def_readwrite("muS", &PyNegfParams::muS)
        .def_readwrite("muD", &PyNegfParams::muD)
        .def_readwrite("isOrthogonal", &PyNegfParams::isOrthogonal)
        .def_readwrite("DCache", &PyNegfParams::DCache)
        .def_readwrite("TCache", &PyNegfParams::TCache)
        .def_readwrite("grcCache", &PyNegfParams::grcCache)
        .def_readwrite("glcCache", &PyNegfParams::glcCache)
        .def_readwrite("GiiCache", &PyNegfParams::GiiCache)
        .def_readwrite("GlCache", &PyNegfParams::GlCache)
        .def_readwrite("GuCache", &PyNegfParams::GuCache)
        .def_readwrite("Gi1Cache", &PyNegfParams::Gi1Cache)
        .def_readwrite("GiNCache", &PyNegfParams::GiNCache)
        .def_readwrite("Giip1Cache", &PyNegfParams::Giip1Cache)
        .def_readwrite("Giim1Cache", &PyNegfParams::Giim1Cache)
        .def(self_ns::str(self_ns::self))
    ;
    
    class_<PyNegfEloop, shared_ptr<PyNegfEloop>, noncopyable>("NegfEloop", 
            init<PyVecGrid&, const PyNegfParams&, const communicator&>())
        .def("run", &PyNegfEloop::run)
        .def("prepare", &PyNegfEloop::prepare)
        .def("preCompute", &PyNegfEloop::preCompute)
        .def("compute", &PyNegfEloop::compute)
        .def("postCompute", &PyNegfEloop::postCompute)
        .def("collect", &PyNegfEloop::collect)
        .def("stepCompleted", &PyNegfEloop::stepCompleted)
        .def("computeTE", &PyNegfEloop::computeTE)
        .def("collectTE", &PyNegfEloop::collectTE)
        .def("enableTE", &PyNegfEloop::enableTE, PyNegfEloop_enableTE())
        .def("saveTE", &PyNegfEloop::saveTE)        
    ;
    
}

}
}