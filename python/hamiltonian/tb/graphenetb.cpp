/* 
 * File:   graphenetb.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 1:04 AM
 * 
 * Tight binding parameters for graphene.
 */

#include "hamiltonian/tb/graphenetb.h"
#include "python/boostpython.hpp"


/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace hamiltonian;

/**
 * Graphene TB parameters.
 */
void export_GrapheneTbParams(){
    double (GrapheneTbParams::*GrapheneTbParams_getec)() = &GrapheneTbParams::ec;
    void (GrapheneTbParams::*GrapheneTbParams_setec)(double) = &GrapheneTbParams::ec;
    double (GrapheneTbParams::*GrapheneTbParams_getdi0)() = &GrapheneTbParams::di0;
    void (GrapheneTbParams::*GrapheneTbParams_setdi0)(double) = &GrapheneTbParams::di0;
    double (GrapheneTbParams::*GrapheneTbParams_getti0)() = &GrapheneTbParams::ti0;
    void (GrapheneTbParams::*GrapheneTbParams_setti0)(double) = &GrapheneTbParams::ti0;
    double (GrapheneTbParams::*GrapheneTbParams_getdo0)() = &GrapheneTbParams::do0;
    void (GrapheneTbParams::*GrapheneTbParams_setdo0)(double) = &GrapheneTbParams::do0;
    double (GrapheneTbParams::*GrapheneTbParams_getto0)() = &GrapheneTbParams::to0;
    void (GrapheneTbParams::*GrapheneTbParams_setto0)(double) = &GrapheneTbParams::to0;
    double (GrapheneTbParams::*GrapheneTbParams_getlmdz)() = &GrapheneTbParams::lmdz;
    void (GrapheneTbParams::*GrapheneTbParams_setlmdz)(double) = &GrapheneTbParams::lmdz;
    double (GrapheneTbParams::*GrapheneTbParams_getlmdxy)() = &GrapheneTbParams::lmdxy;
    void (GrapheneTbParams::*GrapheneTbParams_setlmdxy)(double) = &GrapheneTbParams::lmdxy;
    double (GrapheneTbParams::*GrapheneTbParams_getalpha)() = &GrapheneTbParams::alpha;
    void (GrapheneTbParams::*GrapheneTbParams_setalpha)(double) = &GrapheneTbParams::alpha;
    double (GrapheneTbParams::*GrapheneTbParams_getdoX)() = &GrapheneTbParams::doX;
    void (GrapheneTbParams::*GrapheneTbParams_setdoX)(double) = &GrapheneTbParams::doX;

    class_<GrapheneTbParams, bases<cxhamparams>, shared_ptr<GrapheneTbParams> >(
        "GrapheneTbParams", init<optional<const string &> >())
        .enable_pickling()
        .add_property("ec", GrapheneTbParams_getec, GrapheneTbParams_setec)
        .add_property("di0", GrapheneTbParams_getdi0, GrapheneTbParams_setdi0)
        .add_property("ti0", GrapheneTbParams_getti0, GrapheneTbParams_setti0)   
        .add_property("do0", GrapheneTbParams_getdo0, GrapheneTbParams_setdo0)   
        .add_property("to0", GrapheneTbParams_getto0, GrapheneTbParams_setto0)
        .add_property("lmdz", GrapheneTbParams_getlmdz, GrapheneTbParams_setlmdz)
        .add_property("lmdxy", GrapheneTbParams_getlmdxy, GrapheneTbParams_setlmdxy)
        .add_property("alpha", GrapheneTbParams_getalpha, GrapheneTbParams_setalpha)    
        .add_property("doX", GrapheneTbParams_getdoX, GrapheneTbParams_setdoX)   
    ;
}


}
}

    



