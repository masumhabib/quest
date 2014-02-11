/* 
 * File:   graphenetb.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on January 25, 2014, 1:04 AM
 * 
 * Tight binding parameters for graphene.
 */

#include "graphenetb.h"

namespace qmicad{
using namespace maths::armadillo;

mat GrapheneTbHam::genTwoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    float d, dx, dy, dz;
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    mat hmat =  zeros<mat>(noi, noj);
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
};

mat GrapheneTbHam::genTwoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)
{
    float d, dx, dy, dz;
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    mat smat =  zeros<mat>(noi, noj);
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
            smat = eye(noi, noj);
        }
    }
    return smat;
};


}

