/* 
 * File:   ti.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#include "tikp.h"
#include "../../python/boostpython.hpp"

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
        if (d <= p->dtol){
            hmat = p->eps;
        // nearest neighbor in x
        }else if(abs(d - p->ax) <= p->dtol && abs(dx - p->ax) <= p->dtol){ 
            if (xi > xj){
                hmat = p->t10x;
            }else{
                hmat = p->t01x;
            }
        //nearest neighbor y
        }else if (abs(d - p->ay) <= p->dtol && abs(dy - p->ay) <= p->dtol){
            if(yi > yj){
                hmat = p->t10y;
            }else{
                hmat = p->t01y;
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
        if (d <= p->dtol){
            smat = p->I;
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
    class_<TISurfKpParams, bases<HamParams>, shared_ptr<TISurfKpParams> >("TISurfKpParams")
        .enable_pickling()
        .def_readwrite("ax", &TISurfKpParams::ax)
        .def_readwrite("ay", &TISurfKpParams::ay)
        .def_readwrite("Rx", &TISurfKpParams::Rx)   
        .def_readwrite("Ry", &TISurfKpParams::Ry)   
        .def_readwrite("A2", &TISurfKpParams::A2)
        .def_readwrite("C", &TISurfKpParams::C)
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




