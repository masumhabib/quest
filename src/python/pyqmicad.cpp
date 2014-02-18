/* 
 * File:   pyqmicad.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 3, 2014, 10:58 PM
 * 
 * Python interface to our wrappers.
 * 
 */

#include <python2.6/Python.h>
#include <boost/python.hpp>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>
#include <boost/python/detail/wrap_python.hpp>

#include "pyqmicad.h"
#include "../include/qmicad.hpp"
#include "../utils/vout.h"




namespace qmicad{
namespace python{

/**
 * Prints welcome message.
 */
char const* greet(){   
    static string msg;
    msg  = "      QMICAD: Quantum Mechanics Inspired Computer Aided Design";
    msg += "\n                             v" + qmicad::version + "\n";
    msg += " ------------------------------------------------------------------";
    return msg.c_str();
}

/**
 * Sets application verbosity level.
 */

void setVerbosity(int verb){
    stds::vout.appVerbosity(verbosity(verb));
}

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(LinearPot_VLR, VLR, 3, 5)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableTE, enableTE, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableI1, enableI1, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableI1sx, enableI1, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableI1sy, enableI1, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableI1sz, enableI1, 0, 1)
BOOST_PYTHON_MODULE(qmicad)
{
    using namespace boost::python;
    using boost::mpi::communicator;
    using boost::shared_ptr;
    using boost::noncopyable;
    using namespace utils::stds;
    using namespace maths::geometry;
    
    def("greet", greet);
    def("setVerbosity", setVerbosity);
   
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

    class_<Printable, shared_ptr<Printable> >("Printable", 
            init<optional<const string&> >())
        .def(self_ns::str(self_ns::self))
    ;
    
    /**
     * Spatial vector/position vector.
     */
    class_<PyVec, shared_ptr<PyVec> >("svec", 
            init<optional<double, double, double> >())
        .add_property("X", &PyVec::getx, &PyVec::setx)
        .add_property("Y", &PyVec::gety, &PyVec::sety)
        .add_property("Z", &PyVec::getz, &PyVec::setz)
    ;
    
    /**
     * MPI communicator wrapper.
     */
    class_<Workers, shared_ptr<Workers>, noncopyable>("Workers", 
            init<const communicator&>())
        .def("MyId", &Workers::MyId)
        .def("MasterId", &Workers::MasterId)
        .def("N", &Workers::N)
        .def("AmIMaster", &Workers::AmIMaster)
        .def("IAmMaster", &Workers::IAmMaster)
    ; 
    
    /**
     * Geometry classes
     */
    class_<point, shared_ptr<point> >("Point", 
            init<const double&, const double&>())
        .def_pickle(PointPickler())
    ;
    
    class_<quadrilateral, bases<Printable>, shared_ptr<quadrilateral> >("Quadrilateral", 
            init<const point&, const point&, const point&, 
            const point&, optional<const string&> >())
        .def_pickle(QuadrilateralPickler())
    ;

    /**
     * Wall clock
     */
    class_<Timer, bases<Printable>, shared_ptr<Timer> >("Timer", 
            init<optional<const string&> >())
        .def("tic", &Timer::tic)
        .def("toc", &Timer::toc)
    ;

    
    /**
     * Atom.
     */
    class_<Atom, shared_ptr<Atom> >("Atom")
        .def_pickle(AtomPickler())
    ;
    
    /**
     * Periodic table.
     */
    void (PeriodicTable::*PeriodicTable_add1)(uint, string, uint, uint) = &PeriodicTable::add;
    void (PeriodicTable::*PeriodicTable_add2)(const Atom&) = &PeriodicTable::add;
    class_<PeriodicTable, shared_ptr<PeriodicTable> >("PeriodicTable")
        .def("add", PeriodicTable_add1)
        .def("add", PeriodicTable_add2)
        .def_pickle(PeriodicTablePickler())
    ;

    /**
     * Atomistic geometry of the device.
     */
    class_<AtomicStruct, bases<Printable>, shared_ptr<AtomicStruct> >("AtomicStruct", 
            init<const string&>())
        .def(init<const string&, const PeriodicTable>())
        .def(init<uint, uint, double, double, const PeriodicTable>())
        .def("span", &AtomicStruct::span)
        .def("xmax", &AtomicStruct::xmax)
        .def("xmin", &AtomicStruct::xmin)
        .def("ymax", &AtomicStruct::ymax)
        .def("ymin", &AtomicStruct::ymin)    
        .def("zmax", &AtomicStruct::zmax)
        .def("zmin", &AtomicStruct::zmin)    
    ;
    
    /**
     * Hamiltonian parameters.
     */
    class_<HamParams, bases<Printable>, shared_ptr<HamParams> >("HamParams", 
            init<const string&>())
    ;
    
    /**
     * Graphene k.p parameters.
     */
    class_<GrapheneKpParams, bases<HamParams>, shared_ptr<GrapheneKpParams> >("GrapheneKpParams")
        .def_readwrite("dtol", &GrapheneKpParams::dtol)
        .def_readwrite("ax", &GrapheneKpParams::ax)
        .def_readwrite("ay", &GrapheneKpParams::ay)
        .def_readwrite("K", &GrapheneKpParams::K)   
        .def_readwrite("gamma", &GrapheneKpParams::gamma)
        .def("update", &GrapheneKpParams::update)
    ;

    /**
     * Graphene Hamiltonian.
     */
    void (GrapheneKpHam::*GrapheneKpHam_generate1)(const AtomicStruct&, 
            const AtomicStruct&, uint, uint) = &GrapheneKpHam::generate;
    class_<GrapheneKpHam, shared_ptr<GrapheneKpHam> >("GrapheneKpHam",
            init<const GrapheneKpParams& >())
        .def("setSize", &GrapheneKpHam::setSize)
        .def("setSizeForNegf", &GrapheneKpHam::setSizeForNegf)
        .def("setSizeForBand", &GrapheneKpHam::setSizeForBand)
        .def("H0", &GrapheneKpHam::H0)
        .def("Hl", &GrapheneKpHam::Hl)
        .def("H", &GrapheneKpHam::H)  
        .def("Sl", &GrapheneKpHam::Sl)
        .def("S0", &GrapheneKpHam::S0)
        .def("S", &GrapheneKpHam::S)  
        .def("genDiagBlock", &GrapheneKpHam::genDiagBlock)
        .def("genLowDiagBlock", &GrapheneKpHam::genLowDiagBlock)
        .def("genNearestNeigh", &GrapheneKpHam::genNearestNeigh)
        .def("generate", GrapheneKpHam_generate1)
    ;
    
    /**
     * TI Surface k.p parameters.
     */
    class_<TISurfKpParams, bases<HamParams>, shared_ptr<TISurfKpParams> >("TISurfKpParams")
        .def_readwrite("dtol", &TISurfKpParams::dtol)
        .def_readwrite("ax", &TISurfKpParams::ax)
        .def_readwrite("ay", &TISurfKpParams::ay)
        .def_readwrite("K", &TISurfKpParams::K)   
        .def_readwrite("A2", &TISurfKpParams::A2)
        .def_readwrite("C", &TISurfKpParams::C)
        .def("update", &TISurfKpParams::update)
    ;

    /**
     * TI Surface Hamiltonian.
     */
    void (TISurfKpHam::*TISurfKpHam_generate1)(const AtomicStruct&, 
            const AtomicStruct&, uint, uint) = &TISurfKpHam::generate;    
    class_<TISurfKpHam, shared_ptr<TISurfKpHam> >("TISurfKpHam",
            init<const TISurfKpParams& >())
        .def("setSize", &TISurfKpHam::setSize)
        .def("setSizeForNegf", &TISurfKpHam::setSizeForNegf)
        .def("setSizeForBand", &TISurfKpHam::setSizeForBand)
        .def("H0", &TISurfKpHam::H0)
        .def("Hl", &TISurfKpHam::Hl)
        .def("H", &TISurfKpHam::H)  
        .def("Sl", &TISurfKpHam::Sl)
        .def("S0", &TISurfKpHam::S0)
        .def("S", &TISurfKpHam::S)  
        .def("genDiagBlock", &TISurfKpHam::genDiagBlock)
        .def("genLowDiagBlock", &TISurfKpHam::genLowDiagBlock)
        .def("genNearestNeigh", &TISurfKpHam::genNearestNeigh)
        .def("generate", TISurfKpHam_generate1)
    ;
    
    /**
     * Linear potential
     */  
    shared_ptr<vec> (Potential::*Potential_toOrbPot)(uint, uint) = &Potential::toOrbPot;    
    class_<Potential, bases<Printable>, shared_ptr<Potential> >("Potential", 
            init<const AtomicStruct&, optional<const string&> >())
        .def("addSource", &Potential::addSource)
        .def("addDrain", &Potential::addDrain)
        .def("addGate", &Potential::addGate)
        .def("compute", &Potential::compute)
        .def("exportSvg", &Potential::exportSvg)
        .def("exportPotential", &Potential::exportPotential)
        .def("VD", &Potential::VD)
        .def("VS", &Potential::VS)
        .def("VG", &Potential::VG)
        .def("toOrbPot", Potential_toOrbPot) 
    ;
    
    /**
     * Linear potential
     */  
    class_<LinearPot, bases<Potential>, shared_ptr<LinearPot> >("LinearPot", 
            init<const AtomicStruct&, optional<const string&> >())
        .def("addLinearRegion", &LinearPot::addLinearRegion)
        .def("VLR", &LinearPot::VLR, LinearPot_VLR()) 
    ;
     
    /**
     * Vector grid
     */
    class_<VecGrid, bases<Printable>, shared_ptr<VecGrid> >("VecGrid", 
            init<double, double, double, optional<const string&> >())
        .def(init<optional<double, double, int, const string&> >())
        .def("V", &VecGrid::V)        
        .def("min", &VecGrid::min)
        .def("max", &VecGrid::max)
        .def("del", &VecGrid::del)
        .def("N", &VecGrid::N)
    ;

    /**
     * NEGF parameters.
     */
    class_<NegfParams, bases<Printable>, shared_ptr<NegfParams> >("NegfParams", 
            init<uint>())
        .def("H0", &NegfParams::setH0)
        .def("S0", &NegfParams::setS0)
        .def("Hl", &NegfParams::setHl)
        .def("Sl", &NegfParams::setSl)
        .def("V", &NegfParams::setV)
        .def_readwrite("kT", &NegfParams::kT)
        .def_readwrite("ieta", &NegfParams::ieta)
        .def_readwrite("muS", &NegfParams::muS)
        .def_readwrite("muD", &NegfParams::muD)
        .def_readwrite("isOrthogonal", &NegfParams::isOrthogonal)
        .def_readwrite("DCache", &NegfParams::DCache)
        .def_readwrite("TCache", &NegfParams::TCache)
        .def_readwrite("grcCache", &NegfParams::grcCache)
        .def_readwrite("glcCache", &NegfParams::glcCache)
        .def_readwrite("GiiCache", &NegfParams::GiiCache)
        .def_readwrite("GlCache", &NegfParams::GlCache)
        .def_readwrite("GuCache", &NegfParams::GuCache)
        .def_readwrite("Gi1Cache", &NegfParams::Gi1Cache)
        .def_readwrite("GiNCache", &NegfParams::GiNCache)
        .def_readwrite("Giip1Cache", &NegfParams::Giip1Cache)
        .def_readwrite("Giim1Cache", &NegfParams::Giim1Cache)
    ;
    
    class_<NegfEloop, shared_ptr<NegfEloop>, noncopyable>("NegfEloop", 
            init<VecGrid&, const NegfParams&, const Workers&, 
            optional<bool> >())
        .def("run", &NegfEloop::run)
        .def("save", &NegfEloop::save)
        .def("enableTE", &NegfEloop::enableTE, NegfEloop_enableTE())
        .def("enableI1", &NegfEloop::enableTE, NegfEloop_enableI1())
        .def("enableI1sx", &NegfEloop::enableTE, NegfEloop_enableI1sx())
        .def("enableI1sy", &NegfEloop::enableTE, NegfEloop_enableI1sy())
        .def("enableI1sz", &NegfEloop::enableTE, NegfEloop_enableI1sz())
    ;

//    class_<BandStruct, shared_ptr<BandStruct>, noncopyable>("BandStruct",
//            init<shared_ptr<mat>, const BandStructParams, 
//            const Workers&>()) 
//    ;
    
}

}
}