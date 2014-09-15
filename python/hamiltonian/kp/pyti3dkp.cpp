
/* 
 * File:   pyti3dkp.cpp
 * Author: K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on September 12, 2014, 9:49 AM
 */

#include "TI3DKpParams.h"

namespace qmicad{namespace python{
using namespace hamiltonian;

/**
 * TI Surface k.p parameters.
 */
void export_TI3DKpParams(){
    double (export_TI3DKpParams::*TI3DKpParams_geta)() const = &TI3DKpParams::a;
    void (export_TI3DKpParams::*TI3DKpParams_seta)(double) = &TI3DKpParams::a;
    double (export_TI3DKpParams::*TI3DKpParams_getA1)() const = &TI3DKpParams::A1;
    void (export_TI3DKpParams::*TI3DKpParams_setA1)(double) = &TI3DKpParams::A1;
    double (export_TI3DKpParams::*TI3DKpParams_getA2)() const = &TI3DKpParams::A2;
    void (export_TI3DKpParams::*TI3DKpParams_setA2)(double) = &TI3DKpParams::A2;
    double (export_TI3DKpParams::*TI3DKpParams_getB1)() const = &TI3DKpParams::B1;
    void (export_TI3DKpParams::*TI3DKpParams_setB1)(double) = &TI3DKpParams::B1;
    double (export_TI3DKpParams::*TI3DKpParams_getB2)() const = &TI3DKpParams::B2;
    void (export_TI3DKpParams::*TI3DKpParams_setB2)(double) = &TI3DKpParams::B2;
    double (export_TI3DKpParams::*TI3DKpParams_getC)() const = &TI3DKpParams::C;
    void (export_TI3DKpParams::*TI3DKpParams_setC)(double) = &TI3DKpParams::C;
    double (export_TI3DKpParams::*TI3DKpParams_getD1)() const = &TI3DKpParams::D1;
    void (export_TI3DKpParams::*TI3DKpParams_setD1)(double) = &TI3DKpParams::D1;
    double (export_TI3DKpParams::*TI3DKpParams_getD2)() const = &TI3DKpParams::D2;
    void (export_TI3DKpParams::*TI3DKpParams_setD2)(double) = &TI3DKpParams::D2;
    double (export_TI3DKpParams::*TI3DKpParams_getM)() const = &TI3DKpParams::M;
    void (export_TI3DKpParams::*TI3DKpParams_setM)(double) = &TI3DKpParams::M;

    class_<TI3DKpParams, bases<cxhamparams >, shared_ptr<TI3DKpParams> >("TI3DKpParams",
        init<optional<const string &> >())
        .enable_pickling()
        .add_property("a", TI3DKpParams_geta, TI3DKpParams_seta)
        .add_property("A1", TI3DKpParams_getA1, TI3DKpParams_setA1)
        .add_property("A2", TI3DKpParams_getA2, TI3DKpParams_setA2)
        .add_property("B1", TI3DKpParams_getB1, TI3DKpParams_setB1)
        .add_property("B2", TI3DKpParams_getB2, TI3DKpParams_setB2)
        .add_property("C", TI3DKpParams_getC, TI3DKpParams_setC)
        .add_property("D1", TI3DKpParams_getD1, TI3DKpParams_setD1)
        .add_property("D2", TI3DKpParams_getD2, TI3DKpParams_setD2)
        .add_property("M", TI3DKpParams_getM, TI3DKpParams_setM)
    ;
}


}}




