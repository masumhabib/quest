/* 
 * File:   graphenetb.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on January 25, 2014, 1:04 AM
 * 
 * Tight binding parameters for graphene.
 */

#include "graphenetb.h"
#include "../../python/boostpython.hpp"

namespace qmicad{
namespace hamiltonian{

using namespace maths::armadillo;

cxmat GrapheneTbHam::genTwoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    double d, dx, dy, dz;
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat hmat =  zeros<cxmat>(noi, noj);
    shared_ptr<GrapheneTbParams> p = static_pointer_cast<GrapheneTbParams>(mhp);

    // calculate distance between atom i and atom j
    dx = atomi.X(0) - atomj.X(0);
    dy = atomi.Y(0) - atomj.Y(0);           
    dz = atomi.Z(0) - atomj.Z(0);
    d = sqrt(dx*dx + dy*dy + dz*dz);

    // Assign the the matrix elements based on the distance between 
    // the atoms
    // C-C
    if (atomi.Symbol(0) == "C" && atomj.Symbol(0) == "C"){

        dz = abs(dz); // inter-plane distance

        // site energy
        if (d <= p->dtol){ 
            hmat(0,0) = p->ec;
        // in-plane first nearest neighbor
        }else if (dz <= p->dtol && abs(d - p->di0) <= p->dtol){
            hmat(0,0) = -p->ti0;
        /*// out-of-plane first nearest neighbors
        }else if (abs(delz - do0cc) <= dtol && abs(d - do0cc) <= dtol){
            hmat(ia,ja) = -to0cc;
        }*/
        // out-of-plane some nearest neighbors
        }else if (abs(dz - p->do0) <= p->dtol && abs(d - p->do0) <= (p->doX*p->di0 + p->dtol)){
            // PRB 84, 195421 (2011)
            // DBG hmat(ia,ja) = -to0cc*exp(-3.0*(d - do0cc));

            // PRL 109, 236604 (2012)
            double dxy = sqrt(dx*dx + dy*dy);
            hmat(0,0) = -p->to0*exp(-(d - p->do0)/p->lmdz)*exp(-pow(dxy/p->lmdxy, p->alpha));
        }
    }
    
    return hmat;
};

cxmat GrapheneTbHam::genTwoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    float d, dx, dy, dz;
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat smat =  zeros<cxmat>(noi, noj);
    shared_ptr<GrapheneTbParams> p = static_pointer_cast<GrapheneTbParams>(mhp);

    // calculate distance between atom i and atom j
    dx = atomi.X(0) - atomj.X(0);
    dy = atomi.Y(0) - atomj.Y(0);           
    dz = atomi.Z(0) - atomj.Z(0);
    d = sqrt(dx*dx + dy*dy + dz*dz);

    // Assign the the matrix elements based on the distance between 
    // the atoms
    // C-C
    if (atomi.Symbol(0) == "C" && atomj.Symbol(0) == "C"){
        dz = abs(dz); // inter-plane distance
        // site energy
        if (d <= p->dtol){ 
            smat = eye<cxmat>(noi, noj);
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
 * Graphene TB parameters.
 */
void export_GrapheneTbParams(){
    class_<GrapheneTbParams, bases<HamParams>, shared_ptr<GrapheneTbParams> >("GrapheneTbParams")
        .enable_pickling()
        .def_readwrite("ec", &GrapheneTbParams::ec)
        .def_readwrite("di0", &GrapheneTbParams::di0)
        .def_readwrite("ti0", &GrapheneTbParams::ti0)   
        .def_readwrite("do0", &GrapheneTbParams::do0)   
        .def_readwrite("to0", &GrapheneTbParams::to0)
        .def_readwrite("lmdz", &GrapheneTbParams::lmdz)
        .def_readwrite("lmdxy", &GrapheneTbParams::lmdxy)
        .def_readwrite("alpha", &GrapheneTbParams::alpha)    
        .def_readwrite("doX", &GrapheneTbParams::doX)   
    ;
}

/**
 * Graphene TB Hamiltonian.
 */
void export_GrapheneTbHam(){
    class_<GrapheneTbHam, bases<cxham>, shared_ptr<GrapheneTbHam> >("GrapheneTbHam",
            init<const GrapheneTbParams& >())
    ;
}

}
}

    



