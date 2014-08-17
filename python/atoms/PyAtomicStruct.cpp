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
namespace qmicad{
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
void (PeriodicTable::*PeriodicTable_add1)(uint, string, uint, uint) = &PeriodicTable::add;
void (PeriodicTable::*PeriodicTable_add2)(const Atom&) = &PeriodicTable::add;
void export_PeriodicTable(){
    class_<PeriodicTable, shared_ptr<PeriodicTable> >("PeriodicTable")
        .def("add", PeriodicTable_add1)
        .def("add", PeriodicTable_add2)
        .def_pickle(PeriodicTablePickler())
    ;
}

/**
 * Atomistic geometry of the device.
 */

lvec (AtomicStruct::*AtomicStruct_LatticeVector1)() const = &AtomicStruct::LatticeVector;
void export_AtomicStruct(){
    class_<AtomicStruct, bases<Printable>, shared_ptr<AtomicStruct> >("AtomicStruct", 
            init<>())
        .def(init<const string&>())
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
        .def(self + self)
        .def(self + svec())
        .def(self - svec())
        .def(self + lcoord())
        .def(self - lcoord())
    ;
}

}
}

