/* 
 * File:   ti.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#include "boostpython.hpp"
#include "hamiltonian/kp/tikp.h"


/**
 * Python exporters.
 */
namespace qmicad{namespace python{
using namespace hamiltonian;

    /**
     * TI Surface k.p parameters.
     */
void export_TISurfKpParams(){
    double (TISurfKpParams::*TISurfKpParams_geta)() const = &TISurfKpParams::a;
    void (TISurfKpParams::*TISurfKpParams_seta)(double) = &TISurfKpParams::a;
    double (TISurfKpParams::*TISurfKpParams_getC)() const = &TISurfKpParams::C;
    void (TISurfKpParams::*TISurfKpParams_setC)(double) = &TISurfKpParams::C;
    double (TISurfKpParams::*TISurfKpParams_getA2)() const = &TISurfKpParams::A2;
    void (TISurfKpParams::*TISurfKpParams_setA2)(double) = &TISurfKpParams::A2;
    double (TISurfKpParams::*TISurfKpParams_getK)() const = &TISurfKpParams::K;
    void (TISurfKpParams::*TISurfKpParams_setK)(double) = &TISurfKpParams::K;

    class_<TISurfKpParams, bases<cxhamparams>, shared_ptr<TISurfKpParams> >(
        "TISurfKpParams", init<optional<const string &> >())
        .enable_pickling()
        .add_property("a", TISurfKpParams_geta, TISurfKpParams_seta)
        .add_property("K", TISurfKpParams_getK, TISurfKpParams_setK)   
        .add_property("A2", TISurfKpParams_getA2, TISurfKpParams_setA2)
        .add_property("C", TISurfKpParams_getC, TISurfKpParams_setC)
    ;
}

}}




