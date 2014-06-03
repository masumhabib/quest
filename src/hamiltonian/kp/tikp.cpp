/* 
 * File:   ti.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#include "hamiltonian/kp/tikp.h"
#include "python/boostpython.hpp"

namespace qmicad{
namespace hamiltonian{

cxmat TISurfKpHam::genTwoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
    }

    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat hmat =  zeros<cxmat>(noi, noj);
    shared_ptr<TISurfKpParams> p = static_pointer_cast<TISurfKpParams>(mhp);

    // calculate distance between atom i and atom j
    double xi = atomi.X(0);
    double yi = atomi.Y(0);
    double xj = atomj.X(0);
    double yj = atomj.Y(0);

    double dx = abs(xi - xj);
    double dy = abs(yi - yj);           
    double d = sqrt(dx*dx + dy*dy);

    // Assign the the matrix elements based on the distance between 
    // the lattice points.
    if (atomi.Symbol(0) == "D" && atomj.Symbol(0) == "D"){
        // site energy
        if (d <= p->mdtol){
            hmat = p->meps;
        // nearest neighbor in x
        }else if(abs(d - p->ma) <= p->mdtol && abs(dx - p->ma) <= p->mdtol){ 
            if (xi > xj){
                hmat = p->mt10x;
            }else{
                hmat = p->mt01x;
            }
        //nearest neighbor y
        }else if (abs(d - p->ma) <= p->mdtol && abs(dy - p->ma) <= p->mdtol){
            if(yi > yj){
                hmat = p->mt10y;
            }else{
                hmat = p->mt01y;
            }
        }            
    }
    
    return hmat;
};

cxmat TISurfKpHam::genTwoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
    }

    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat smat =  zeros<cxmat>(noi, noj);
    shared_ptr<TISurfKpParams> p = static_pointer_cast<TISurfKpParams>(mhp);

    // calculate distance between atom i and atom j
    double xi = atomi.X(0);
    double yi = atomi.Y(0);
    double xj = atomj.X(0);
    double yj = atomj.Y(0);

    double dx = abs(xi - xj);
    double dy = abs(yi - yj);           
    double d = sqrt(dx*dx + dy*dy);

    // Assign the the matrix elements based on the distance between 
    // the lattice points.
    if (atomi.Symbol(0) == "D" && atomj.Symbol(0) == "D"){
        // site energy
        if (d <= p->mdtol){
            smat = p->mI;
        }
    }   
    return smat;
};

}
}

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace hamiltonian;

    /**
     * TI Surface k.p parameters.
     */
void export_TISurfKpParams(){
    double (TISurfKpParams::*TISurfKpParams_geta)() = &TISurfKpParams::a;
    void (TISurfKpParams::*TISurfKpParams_seta)(double) = &TISurfKpParams::a;
    double (TISurfKpParams::*TISurfKpParams_getC)() = &TISurfKpParams::C;
    void (TISurfKpParams::*TISurfKpParams_setC)(double) = &TISurfKpParams::C;
    double (TISurfKpParams::*TISurfKpParams_getA2)() = &TISurfKpParams::A2;
    void (TISurfKpParams::*TISurfKpParams_setA2)(double) = &TISurfKpParams::A2;
    double (TISurfKpParams::*TISurfKpParams_getK)() = &TISurfKpParams::K;
    void (TISurfKpParams::*TISurfKpParams_setK)(double) = &TISurfKpParams::K;

    class_<TISurfKpParams, bases<HamParams>, shared_ptr<TISurfKpParams> >("TISurfKpParams")
        .enable_pickling()
        .add_property("a", TISurfKpParams_geta, TISurfKpParams_seta)
        .add_property("K", TISurfKpParams_getK, TISurfKpParams_setK)   
        .add_property("A2", TISurfKpParams_getA2, TISurfKpParams_setA2)
        .add_property("C", TISurfKpParams_getC, TISurfKpParams_setC)
    ;
}

    /**
     * TI Surface Hamiltonian.
     */
void export_TISurfKpHam(){
    class_<TISurfKpHam, bases<cxham>, shared_ptr<TISurfKpHam> >("TISurfKpHam",
            init<const TISurfKpParams& >())
    ;
}

}
}




