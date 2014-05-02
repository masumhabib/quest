/* 
 * File:   graphenekp.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */


#include "graphenekp.h"
#include "../../python/boostpython.hpp"

namespace qmicad{
namespace hamiltonian{


cxmat GrapheneKpHam::genTwoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj)
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1) {
        throw invalid_argument("GrapheneHamGen(): atomi and atomj must should contain one atom each.");
    }
    
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat hmat =  zeros<cxmat>(noi, noj);
    shared_ptr<GrapheneKpParams> p = static_pointer_cast<GrapheneKpParams>(mhp);

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

cxmat GrapheneKpHam::genTwoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj)
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1) {
        throw invalid_argument("GrapheneHamGen(): atomi and atomj must should contain one atom each.");
    }
    
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat smat =  zeros<cxmat>(noi, noj);
    shared_ptr<GrapheneKpParams> p = static_pointer_cast<GrapheneKpParams>(mhp);

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
        // overlap matrix of atom i.
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
 * Graphene k.p parameters.
 */
void export_GrapheneKpParams(){
    double (GrapheneKpParams::*GrapheneKpParams_geta)() = &GrapheneKpParams::a;
    void (GrapheneKpParams::*GrapheneKpParams_seta)(double) = &GrapheneKpParams::a;
    double (GrapheneKpParams::*GrapheneKpParams_getgamma)() = &GrapheneKpParams::gamma;
    void (GrapheneKpParams::*GrapheneKpParams_setgamma)(double) = &GrapheneKpParams::gamma;
    double (GrapheneKpParams::*GrapheneKpParams_getK)() = &GrapheneKpParams::K;
    void (GrapheneKpParams::*GrapheneKpParams_setK)(double) = &GrapheneKpParams::K;    
    class_<GrapheneKpParams, bases<HamParams>, shared_ptr<GrapheneKpParams> >("GrapheneKpParams")
        .enable_pickling()
        .add_property("a", GrapheneKpParams_geta, GrapheneKpParams_seta)
        .add_property("K", GrapheneKpParams_getK, GrapheneKpParams_setK)   
        .add_property("gamma", GrapheneKpParams_getgamma, GrapheneKpParams_setgamma)
    ;
}

/**
 * Graphene k.p Hamiltonian.
 */
void export_GrapheneKpHam(){
    class_<GrapheneKpHam, bases<cxham >, shared_ptr<GrapheneKpHam> >("GrapheneKpHam",
            init<const GrapheneKpParams& >())
    ;
}

}
}



