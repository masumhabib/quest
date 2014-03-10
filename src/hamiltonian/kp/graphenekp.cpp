/* 
 * File:   graphenekp.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
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
 * Graphene k.p parameters.
 */
void export_GrapheneKpParams(){
    class_<GrapheneKpParams, bases<HamParams>, shared_ptr<GrapheneKpParams> >("GrapheneKpParams")
        .def_readwrite("ax", &GrapheneKpParams::ax)
        .def_readwrite("ay", &GrapheneKpParams::ay)
        .def_readwrite("Rx", &GrapheneKpParams::Rx)   
        .def_readwrite("Ry", &GrapheneKpParams::Ry)   
        .def_readwrite("gamma", &GrapheneKpParams::gamma)
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



