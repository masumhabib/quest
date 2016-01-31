/* 
 * File:   graphenetb.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 1:04 AM
 * 
 * Tight binding parameters for graphene.
 */

#include "hamiltonian/tb/graphenetb.h"

namespace qmicad{
namespace hamiltonian{

using namespace maths::armadillo;

// Constructor
GrapheneTbParams::GrapheneTbParams(string prefix):
    cxhamparams(prefix)
{
    mTitle = "Graphene tight binding parameters";
    setDefaultParams();
}
        
string GrapheneTbParams::toString() const { 
    stringstream ss;
    ss << cxhamparams::toString() << ": " << endl;
    ss << mPrefix << " dtol  = " << mdtol << endl;
    ss << mPrefix << " acc   = " << mdi0 << endl;
    ss << mPrefix << " ec    = " << mec << endl;
    ss << mPrefix << " di0   = " << mdi0 << endl;
    ss << mPrefix << " ti0   = " << mti0 << endl;
    ss << mPrefix << " do0   = " << mdo0 << endl;
    ss << mPrefix << " to0   = " << mto0 << endl;
    ss << mPrefix << " doX   = " << mdoX << endl;
    ss << mPrefix << " lmdz  = " << mlmdz << endl;
    ss << mPrefix << " lmdxy = " << mlmdxy << endl; 
    ss << mPrefix << " alpha = " << malpha << endl;

    return ss.str(); 
};   

void GrapheneTbParams::update(){
}

cxmat GrapheneTbParams::twoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj) const
{
    double d, dx, dy, dz;
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat hmat =  zeros<cxmat>(noi, noj);

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
        if (d <= mdtol){ 
            hmat(0,0) = mec;
        // in-plane first nearest neighbor
        }else if (dz <= mdtol && abs(d - mdi0) <= mdtol){
            hmat(0,0) = -mti0;
        /*// out-of-plane first nearest neighbors
        }else if (abs(delz - do0cc) <= dtol && abs(d - do0cc) <= dtol){
            hmat(ia,ja) = -to0cc;
        }*/
        // out-of-plane some nearest neighbors
        }else if (abs(dz - mdo0) <= mdtol && abs(d - mdo0) <= (mdoX*mdi0 + mdtol)){
            // PRB 84, 195421 (2011)
            // DBG hmat(ia,ja) = -to0cc*exp(-3.0*(d - do0cc));

            // PRL 109, 236604 (2012)
            double dxy = sqrt(dx*dx + dy*dy);
            hmat(0,0) = -mto0*exp(-(d - mdo0)/mlmdz)*exp(-pow(dxy/mlmdxy, malpha));
        }
    }
    
    return hmat;
};

cxmat GrapheneTbParams::twoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj) const
{
    float d, dx, dy, dz;
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat smat =  zeros<cxmat>(noi, noj);

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
        if (d <= mdtol){ 
            smat = eye<cxmat>(noi, noj);
        }
    }
    
    return smat;
};


}
}

    



