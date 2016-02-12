/* 
 * File:   graphenekp.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */


#include "hamiltonian/kp/graphenekp.h"
#include "boostpython.hpp"

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace hamiltonian;

/**
 * Graphene k.p parameters.
 */
void export_GrapheneKpParams(){

    class_<DiracKpParams, bases<cxhamparams>, shared_ptr<DiracKpParams> >(
        "DiracKpParams", init<optional<const string &> >())
    ;

    double (GrapheneKpParams::*GrapheneKpParams_geta)() const = &GrapheneKpParams::a;
    void (GrapheneKpParams::*GrapheneKpParams_seta)(double) = &GrapheneKpParams::a;
    double (GrapheneKpParams::*GrapheneKpParams_getgamma)() const = &GrapheneKpParams::gamma;
    void (GrapheneKpParams::*GrapheneKpParams_setgamma)(double) = &GrapheneKpParams::gamma;
    double (GrapheneKpParams::*GrapheneKpParams_getK)() const = &GrapheneKpParams::K;
    void (GrapheneKpParams::*GrapheneKpParams_setK)(double) = &GrapheneKpParams::K;    
    class_<GrapheneKpParams, bases<DiracKpParams>, shared_ptr<GrapheneKpParams> >(
        "GrapheneKpParams", init<optional<const string &> >())
        .enable_pickling()
        .add_property("a", GrapheneKpParams_geta, GrapheneKpParams_seta)
        .add_property("K", GrapheneKpParams_getK, GrapheneKpParams_setK)   
        .add_property("gamma", GrapheneKpParams_getgamma, GrapheneKpParams_setgamma)
    ;

    class_<GrapheneOneValleyKpParams, bases<GrapheneKpParams>, shared_ptr<GrapheneOneValleyKpParams> >(
        "GrapheneOneValleyKpParams", init<optional<const string &> >())
        .enable_pickling()
    ;

    class_<GrapheneTwoValleyKpParams, bases<GrapheneKpParams>, shared_ptr<GrapheneTwoValleyKpParams> >(
        "GrapheneTwoValleyKpParams", init<optional<const string &> >())
        .enable_pickling()
    ;

}

}
}



