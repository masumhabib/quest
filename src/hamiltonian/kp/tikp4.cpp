/* 
 * File:   tikp4.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 28, 2014, 08:31 PM
 */

#include "tikp4.h"
#include "../../python/boostpython.hpp"

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
void export_TISurfKpParams4(){
    class_<TISurfKpParams4, bases<HamParams>, shared_ptr<TISurfKpParams4> >("TISurfKpParams4")
        .enable_pickling()
        .def_readwrite("ax", &TISurfKpParams4::ax)
        .def_readwrite("ay", &TISurfKpParams4::ay)
        .def_readwrite("Rx", &TISurfKpParams4::Rx)   
        .def_readwrite("Ry", &TISurfKpParams4::Ry)   
        .def_readwrite("A2", &TISurfKpParams4::A2)
        .def_readwrite("C", &TISurfKpParams4::C)
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




