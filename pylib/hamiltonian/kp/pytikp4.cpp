/* 
 * File:   tikp4.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 28, 2014, 08:31 PM
 */

#include "boostpython.hpp"
#include "hamiltonian/kp/tikp4.h"


/**
 * Python exporters.
 */
namespace qmicad{namespace python{
using namespace hamiltonian;

/**
 * TI Surface k.p parameters.
 */
void export_TISurfKpParams4(){
    double (TISurfKpParams4::*TISurfKpParams4_geta)() const = &TISurfKpParams4::a;
    void (TISurfKpParams4::*TISurfKpParams4_seta)(double) = &TISurfKpParams4::a;
    double (TISurfKpParams4::*TISurfKpParams4_getC)() const = &TISurfKpParams4::C;
    void (TISurfKpParams4::*TISurfKpParams4_setC)(double) = &TISurfKpParams4::C;
    double (TISurfKpParams4::*TISurfKpParams4_getA2)() const = &TISurfKpParams4::A2;
    void (TISurfKpParams4::*TISurfKpParams4_setA2)(double) = &TISurfKpParams4::A2;
    double (TISurfKpParams4::*TISurfKpParams4_getK)() const = &TISurfKpParams4::K;
    void (TISurfKpParams4::*TISurfKpParams4_setK)(double) = &TISurfKpParams4::K;
    class_<TISurfKpParams4, bases<cxhamparams >, shared_ptr<TISurfKpParams4> >("TISurfKpParams4",
        init<optional<const string &> >())
        .enable_pickling()
        .add_property("a", TISurfKpParams4_geta, TISurfKpParams4_seta)
        .add_property("K", TISurfKpParams4_getK, TISurfKpParams4_setK)   
        .add_property("A2", TISurfKpParams4_getA2, TISurfKpParams4_setA2)
        .add_property("C", TISurfKpParams4_getC, TISurfKpParams4_setC)
    ;
}


}}




