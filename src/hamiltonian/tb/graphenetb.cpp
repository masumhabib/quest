/* 
 * File:   graphenetb.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 1:04 AM
 * 
 * Tight binding parameters for graphene.
 */

#include "hamiltonian/tb/graphenetb.h"
#include "python/boostpython.hpp"

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
    double (GrapheneTbParams::*GrapheneTbParams_getec)() = &GrapheneTbParams::ec;
    void (GrapheneTbParams::*GrapheneTbParams_setec)(double) = &GrapheneTbParams::ec;
    double (GrapheneTbParams::*GrapheneTbParams_getdi0)() = &GrapheneTbParams::di0;
    void (GrapheneTbParams::*GrapheneTbParams_setdi0)(double) = &GrapheneTbParams::di0;
    double (GrapheneTbParams::*GrapheneTbParams_getti0)() = &GrapheneTbParams::ti0;
    void (GrapheneTbParams::*GrapheneTbParams_setti0)(double) = &GrapheneTbParams::ti0;
    double (GrapheneTbParams::*GrapheneTbParams_getdo0)() = &GrapheneTbParams::do0;
    void (GrapheneTbParams::*GrapheneTbParams_setdo0)(double) = &GrapheneTbParams::do0;
    double (GrapheneTbParams::*GrapheneTbParams_getto0)() = &GrapheneTbParams::to0;
    void (GrapheneTbParams::*GrapheneTbParams_setto0)(double) = &GrapheneTbParams::to0;
    double (GrapheneTbParams::*GrapheneTbParams_getlmdz)() = &GrapheneTbParams::lmdz;
    void (GrapheneTbParams::*GrapheneTbParams_setlmdz)(double) = &GrapheneTbParams::lmdz;
    double (GrapheneTbParams::*GrapheneTbParams_getlmdxy)() = &GrapheneTbParams::lmdxy;
    void (GrapheneTbParams::*GrapheneTbParams_setlmdxy)(double) = &GrapheneTbParams::lmdxy;
    double (GrapheneTbParams::*GrapheneTbParams_getalpha)() = &GrapheneTbParams::alpha;
    void (GrapheneTbParams::*GrapheneTbParams_setalpha)(double) = &GrapheneTbParams::alpha;
    double (GrapheneTbParams::*GrapheneTbParams_getdoX)() = &GrapheneTbParams::doX;
    void (GrapheneTbParams::*GrapheneTbParams_setdoX)(double) = &GrapheneTbParams::doX;

    class_<GrapheneTbParams, bases<cxhamparams>, shared_ptr<GrapheneTbParams> >(
        "GrapheneTbParams", init<optional<const string &> >())
        .enable_pickling()
        .add_property("ec", GrapheneTbParams_getec, GrapheneTbParams_setec)
        .add_property("di0", GrapheneTbParams_getdi0, GrapheneTbParams_setdi0)
        .add_property("ti0", GrapheneTbParams_getti0, GrapheneTbParams_setti0)   
        .add_property("do0", GrapheneTbParams_getdo0, GrapheneTbParams_setdo0)   
        .add_property("to0", GrapheneTbParams_getto0, GrapheneTbParams_setto0)
        .add_property("lmdz", GrapheneTbParams_getlmdz, GrapheneTbParams_setlmdz)
        .add_property("lmdxy", GrapheneTbParams_getlmdxy, GrapheneTbParams_setlmdxy)
        .add_property("alpha", GrapheneTbParams_getalpha, GrapheneTbParams_setalpha)    
        .add_property("doX", GrapheneTbParams_getdoX, GrapheneTbParams_setdoX)   
    ;
}


}
}

    



