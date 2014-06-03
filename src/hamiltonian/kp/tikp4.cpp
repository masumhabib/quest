/* 
 * File:   tikp4.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 28, 2014, 08:31 PM
 */

#include "hamiltonian/kp/tikp4.h"
#include "python/boostpython.hpp"

namespace qmicad{
namespace hamiltonian{

cxmat TISurfKpHam4::genTwoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
    }

    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat hmat =  zeros<cxmat>(noi, noj);
    shared_ptr<TISurfKpParams4> p = static_pointer_cast<TISurfKpParams4>(mhp);

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

cxmat TISurfKpHam4::genTwoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
    }

    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat smat =  zeros<cxmat>(noi, noj);
    shared_ptr<TISurfKpParams4> p = static_pointer_cast<TISurfKpParams4>(mhp);

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
void export_TISurfKpParams4(){
    double (TISurfKpParams4::*TISurfKpParams4_geta)() = &TISurfKpParams4::a;
    void (TISurfKpParams4::*TISurfKpParams4_seta)(double) = &TISurfKpParams4::a;
    double (TISurfKpParams4::*TISurfKpParams4_getC)() = &TISurfKpParams4::C;
    void (TISurfKpParams4::*TISurfKpParams4_setC)(double) = &TISurfKpParams4::C;
    double (TISurfKpParams4::*TISurfKpParams4_getA2)() = &TISurfKpParams4::A2;
    void (TISurfKpParams4::*TISurfKpParams4_setA2)(double) = &TISurfKpParams4::A2;
    double (TISurfKpParams4::*TISurfKpParams4_getK)() = &TISurfKpParams4::K;
    void (TISurfKpParams4::*TISurfKpParams4_setK)(double) = &TISurfKpParams4::K;
    class_<TISurfKpParams4, bases<HamParams>, shared_ptr<TISurfKpParams4> >("TISurfKpParams4")
        .enable_pickling()
        .add_property("a", TISurfKpParams4_geta, TISurfKpParams4_seta)
        .add_property("K", TISurfKpParams4_getK, TISurfKpParams4_setK)   
        .add_property("A2", TISurfKpParams4_getA2, TISurfKpParams4_setA2)
        .add_property("C", TISurfKpParams4_getC, TISurfKpParams4_setC)
    ;
}

    /**
     * TI Surface Hamiltonian.
     */
void export_TISurfKpHam4(){
    class_<TISurfKpHam4, bases<cxham>, shared_ptr<TISurfKpHam4> >("TISurfKpHam4",
            init<const TISurfKpParams4& >())
    ;
}

}
}




