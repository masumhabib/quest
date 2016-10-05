/* 
 * File:   AtomicStruct.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 5, 2013, 1:27 PM
 * 
 * Description: The class Atoms encapsulates everything about the atoms present in
 * the device structure.
 * 
 */

#include "PyAtomicStruct.h"

/**
 * Python exporters.
 */
namespace quest{
namespace python{
using namespace atoms;
using namespace utils::stds;

/**
 * Atom.
 */
void export_Atom(){
    class_<Atom, shared_ptr<Atom> >("Atom", "Represents an Atom.")
        .def_pickle(AtomPickler())
    ;
}

/**
 * Periodic table.
 */
Atom (PeriodicTable::*PeriodicTable_getitem)(uint) const = &PeriodicTable::operator [];
void (PeriodicTable::*PeriodicTable_add1)(uint, string, uint, uint) = &PeriodicTable::add;
void (PeriodicTable::*PeriodicTable_add2)(const Atom&) = &PeriodicTable::add;
void export_PeriodicTable(){
    class_<PeriodicTable, shared_ptr<PeriodicTable> >("PeriodicTable")
        .def("add", PeriodicTable_add1)
        .def("add", PeriodicTable_add2)
        .def("__getitem__", PeriodicTable_getitem)
        .def_pickle(PeriodicTablePickler())
    ;
}

/**
 * Atomistic geometry of the device.
 */
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(AtomicStruct_genSimpleCubicStruct1, genSimpleCubicStruct, 3, 5)
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(AtomicStruct_genGNR1, genGNR, 4, 5);
void (AtomicStruct::*AtomicStruct_genSimpleCubicStruct1)(const Atom &, double, uint, uint, uint) = &AtomicStruct::genSimpleCubicStruct;
//void (AtomicStruct::*AtomicStruct_genSimpleCubicStruct2)(const Atom &, double, double, double, double) = &AtomicStruct::genSimpleCubicStruct;
void (AtomicStruct::*AtomicStruct_genGNR1)(const Atom &, double, uint, uint, uint) = &AtomicStruct::genGNR;
//void (AtomicStruct::*AtomicStruct_genGNR2)(const Atom &, double, double, double, double) = &AtomicStruct::genGNR;
lvec (AtomicStruct::*AtomicStruct_LatticeVector1)() const = &AtomicStruct::LatticeVector;
void export_AtomicStruct(){
    class_<AtomicStruct, bases<Printable>, shared_ptr<AtomicStruct> >("AtomicStruct", 
            init<>())
        .def(init<const string&>())
        .def(init<const PeriodicTable>())
        .def(init<const string&, const PeriodicTable>())
        .def_pickle(AtomicStructPickler())
        .add_property("xmax", &AtomicStruct::xmax)
        .add_property("xmin", &AtomicStruct::xmin)
        .add_property("xl", &AtomicStruct::xl)
        .add_property("ymax", &AtomicStruct::ymax)
        .add_property("ymin", &AtomicStruct::ymin)    
        .add_property("yl", &AtomicStruct::yl)
        .add_property("zmax", &AtomicStruct::zmax)
        .add_property("zmin", &AtomicStruct::zmin) 
        .add_property("zl", &AtomicStruct::zl)
        .add_property("LatticeVector", AtomicStruct_LatticeVector1) 
        .add_property("NumOfAtoms", &AtomicStruct::NumOfAtoms) 
        .add_property("NumOfOrbitals", &AtomicStruct::NumOfOrbitals) 
        .add_property("NumOfElectrons", &AtomicStruct::NumOfElectrons) 
        .def("span", &AtomicStruct::span)
        .def("genRectLattAtoms", &AtomicStruct::genRectLattAtoms)
        .def("genSimpleCubicStruct", AtomicStruct_genSimpleCubicStruct1)
        //.def("genSimpleCubicStruct", AtomicStruct_genSimpleCubicStruct2)
        .def("genGNR", AtomicStruct_genGNR1)
        //.def("genGNR", AtomicStruct_genGNR2)
        .def("exportGjf", &AtomicStruct::exportGjf)
        .def("importGjf", &AtomicStruct::importGjf)
        .def(self + self)
        .def(self + svec())
        .def(self - svec())
        .def(self + lcoord())
        .def(self - lcoord())
    ;
}

}
}

