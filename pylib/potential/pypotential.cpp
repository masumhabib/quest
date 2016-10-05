/* 
 * File:   potential.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 27, 2014, 3:47 PM
 */

#include "boostpython.hpp"
#include "potential/potential.h"


/**
 * Python exporters.
 */
namespace quest{
namespace python{
using namespace potential;

/**
 * Linear potential
 */  
vec (Potential::*Potential_toOrbPot)(uint, uint) = &Potential::toOrbPot;
double (Potential::*Potential_Vatom1)(uint) = &Potential::Vatom;
void (Potential::*Potential_Vatom2)(uint, double) = &Potential::Vatom;
void export_Potential(){    
    class_<Potential, bases<Printable>, shared_ptr<Potential> >("Potential", 
            init<optional<AtomicStruct::ptr, const string&> >())
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
        .def("Vatom", Potential_Vatom1) 
        .def("Vatom", Potential_Vatom2) 
        .add_property("NG", &Potential::NG)
    ;
}

}
}
    

