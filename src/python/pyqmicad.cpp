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
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableI, enableI, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ham_H, H, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(ham_S, S, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxham_H, H, 1, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(cxham_S, S, 1, 2)
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
    class_<svec, shared_ptr<svec> >("svec")
    ;

    /**
     * Position vector. Just a wrapper of svec.
     */
    class_<PyVec, bases<svec>, shared_ptr<PyVec> >("pvec", 
            init<optional<double, double, double> >())
        .def(init<const svec&>())    
        .add_property("X", &PyVec::getx, &PyVec::setx)
        .add_property("Y", &PyVec::gety, &PyVec::sety)
        .add_property("Z", &PyVec::getz, &PyVec::setz)
    ;
    
    /**
     * Lattice vector.
     */
    class_<lvec, bases<Printable>, shared_ptr<lvec> >("lvec", 
            init<optional<const string&> >())
        .def_readwrite("a1", &lvec::a1)
        .def_readwrite("a2", &lvec::a2)
        .def_readwrite("a3", &lvec::a3)
        .def("la1", &lvec::la1)
        .def("la2", &lvec::la2)
        .def("la3", &lvec::la3)
    ;
    /**
     * Lattice coordinate.
     */
    class_<lcoord, bases<Printable>, shared_ptr<lcoord> >("lcoord", 
            init<int, int, int, optional<const string&> >())
        .def_readwrite("n1", &lcoord::n1)
        .def_readwrite("n2", &lcoord::n2)
        .def_readwrite("n3", &lcoord::n3)
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
    lvec (AtomicStruct::*AtomicStruct_LatticeVector1)() const = &AtomicStruct::LatticeVector;
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
        .def("LatticeVector", AtomicStruct_LatticeVector1) 
        .def("NumOfElectrons", &AtomicStruct::NumOfElectrons) 
        .def("NumOfOrbitals", &AtomicStruct::NumOfOrbitals) 
        .def(self + svec())
        .def(self - svec())
        .def(self + lcoord())
        .def(self - lcoord())
    ;
    
    /**
     * Hamiltonian parameters.
     */
    class_<HamParams, bases<Printable>, shared_ptr<HamParams> >("HamParams", 
            init<const string&>())
        .def_readwrite("dtol", &HamParams::dtol)
        .def("update", &HamParams::update)
    ;

    /**
     * Real Hamiltonian base.
     */
    void (ham::*ham_generate1)(const AtomicStruct&, 
            const AtomicStruct&, uint, uint) = &ham::generate;
    class_<ham, shared_ptr<ham > >("Hamiltonian", no_init)
        .def("setSize", &ham::setSize)
        .def("setSizeForNegf", &ham::setSizeForNegf)
        .def("setSizeForBand", &ham::setSizeForBand)
        .def("H0", &ham::H0)
        .def("Hl", &ham::Hl)
        .def("H", &ham::H, ham_H())  
        .def("Sl", &ham::Sl)
        .def("S0", &ham::S0)
        .def("S", &ham::S, ham_S())  
        .def("genDiagBlock", &ham::genDiagBlock)
        .def("genLowDiagBlock", &ham::genLowDiagBlock)
        .def("genNearestNeigh", &ham::genNearestNeigh)
        .def("generate", ham_generate1)
    ;
    
    /**
     * Complex Hamiltonian base.
     */
    void (cxham::*cxham_generate1)(const AtomicStruct&, 
            const AtomicStruct&, uint, uint) = &cxham::generate;
    class_<cxham, shared_ptr<cxham > >("CxHamiltonian", no_init)
        .def("setSize", &cxham::setSize)
        .def("setSizeForNegf", &cxham::setSizeForNegf)
        .def("setSizeForBand", &cxham::setSizeForBand)
        .def("H0", &cxham::H0)
        .def("Hl", &cxham::Hl)
        .def("H", &cxham::H, cxham_H())  
        .def("Sl", &cxham::Sl)
        .def("S0", &cxham::S0)
        .def("S", &cxham::S, cxham_S())  
        .def("genDiagBlock", &cxham::genDiagBlock)
        .def("genLowDiagBlock", &cxham::genLowDiagBlock)
        .def("genNearestNeigh", &cxham::genNearestNeigh)
        .def("generate", cxham_generate1)
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
        .def_readwrite("DCache", &NegfParams::DCacheEnabled)
        .def_readwrite("TCache", &NegfParams::TCacheEnabled)
        .def_readwrite("grcCache", &NegfParams::grcCacheEnabled)
        .def_readwrite("glcCache", &NegfParams::glcCacheEnabled)
        .def_readwrite("GiiCache", &NegfParams::GiiCacheEnabled)
        .def_readwrite("Gi1Cache", &NegfParams::Gi1CacheEnabled)
        .def_readwrite("GiNCache", &NegfParams::GiNCacheEnabled)
        .def_readwrite("Giip1Cache", &NegfParams::Giip1CacheEnabled)
        .def_readwrite("Giim1Cache", &NegfParams::Giim1CacheEnabled)
    ;
    
    class_<NegfEloop, shared_ptr<NegfEloop>, noncopyable>("NegfEloop", 
            init<VecGrid&, const NegfParams&, const Workers&, 
            optional<bool> >())
        .def("run", &NegfEloop::run)
        .def("save", &NegfEloop::save)
        .def("enableTE", &NegfEloop::enableTE, NegfEloop_enableTE())
        .def("enableI", &NegfEloop::enableI, NegfEloop_enableI())
    ;

    class_<KPoints, bases<Printable>, shared_ptr<KPoints>, noncopyable>("KPoints",
            init<optional<const string&> >()) 
        .def("addKPoint", &KPoints::addKPoint)
        .def("addKLine", &KPoints::addKLine)
        .def("addKRect", &KPoints::addKRect)
        .def("kp", &KPoints::kp)
        .def("N", &KPoints::N)
    ;

    class_<BandStructParams, bases<Printable>, shared_ptr<BandStructParams>, noncopyable>("BandStructParams",
        init<uint, optional<const string&> >()) 
        .def("H", &BandStructParams::setH)
        .def("S", &BandStructParams::setS)
        .def("lc", &BandStructParams::setLattCoord)
        .def_readwrite("nb", &BandStructParams::nb)
        .def_readwrite("ne", &BandStructParams::ne)
        .def_readwrite("no", &BandStructParams::no)
        .def_readwrite("lv", &BandStructParams::lv)
        .def_readwrite("isOrthogonal", &BandStructParams::isOrthogonal)
    ;

    class_<BandStruct, shared_ptr<BandStruct>, noncopyable>("BandStruct",
            init<shared_ptr<mat>, const BandStructParams, 
            const Workers&>())
        .def("run", &BandStruct::run)
        .def("save", &BandStruct::save)
    ;
    
    /**
     * Graphene TB parameters.
     */
    class_<GrapheneTbParams, bases<HamParams>, shared_ptr<GrapheneTbParams> >("GrapheneTbParams")
        .def_readwrite("ec", &GrapheneTbParams::ec)
        .def_readwrite("di0", &GrapheneTbParams::di0)
        .def_readwrite("ti0", &GrapheneTbParams::ti0)   
        .def_readwrite("do0", &GrapheneTbParams::do0)   
        .def_readwrite("to0", &GrapheneTbParams::to0)
        .def_readwrite("lmdz", &GrapheneTbParams::lmdz)
        .def_readwrite("lmdxy", &GrapheneTbParams::lmdxy)
        .def_readwrite("alpha", &GrapheneTbParams::alpha)    
        .def_readwrite("doX", &GrapheneTbParams::doX)   
    ;

    /**
     * Graphene TB Hamiltonian.
     */
    class_<GrapheneTbHam, bases<cxham>, shared_ptr<GrapheneTbHam> >("GrapheneTbHam",
            init<const GrapheneTbParams& >())
    ;
    
    /**
     * Graphene k.p parameters.
     */
    class_<GrapheneKpParams, bases<HamParams>, shared_ptr<GrapheneKpParams> >("GrapheneKpParams")
        .def_readwrite("ax", &GrapheneKpParams::ax)
        .def_readwrite("ay", &GrapheneKpParams::ay)
        .def_readwrite("Rx", &GrapheneKpParams::Rx)   
        .def_readwrite("Ry", &GrapheneKpParams::Ry)   
        .def_readwrite("gamma", &GrapheneKpParams::gamma)
    ;

    /**
     * Graphene k.p Hamiltonian.
     */
    class_<GrapheneKpHam, bases<cxham >, shared_ptr<GrapheneKpHam> >("GrapheneKpHam",
            init<const GrapheneKpParams& >())
    ;
    
    /**
     * TI Surface k.p parameters.
     */
    class_<TISurfKpParams, bases<HamParams>, shared_ptr<TISurfKpParams> >("TISurfKpParams")
        .def_readwrite("ax", &TISurfKpParams::ax)
        .def_readwrite("ay", &TISurfKpParams::ay)
        .def_readwrite("Rx", &TISurfKpParams::Rx)   
        .def_readwrite("Ry", &TISurfKpParams::Ry)   
        .def_readwrite("A2", &TISurfKpParams::A2)
        .def_readwrite("C", &TISurfKpParams::C)
    ;

    /**
     * TI Surface Hamiltonian.
     */
    class_<TISurfKpHam, bases<cxham>, shared_ptr<TISurfKpHam> >("TISurfKpHam",
            init<const TISurfKpParams& >())
    ;
    
    
}

}
}