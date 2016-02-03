/* 
 * File:   graphenekp.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */


#include <maths/svec.h>

#include "hamiltonian/kp/graphenekp.h"

namespace qmicad{
namespace hamiltonian{

GrapheneKpParams::GrapheneKpParams(const string &prefix):cxhamparams(prefix)
{
    mTitle  = "Graphene k.p parameters";
    mI      = eye<cxmat>(2,2);    
    
    setDefaultParams();
    update();
}

string GrapheneKpParams::toString() const { 
    stringstream ss;
    ss << cxhamparams::toString() << ":" << endl;
    ss << mPrefix << " a     = " << ma << endl;
    ss << mPrefix << " K     = " << mK  << endl;
    ss << mPrefix << " gamma = " << mgamma << endl;
    ss << mPrefix << " dtol  = " << mdtol << endl;
    ss << mPrefix << " Bz    = " << mBz << endl;
    ss << mPrefix << " A     = " << (mBzGauge == coord::X?" (-Bz*y, 0, 0)":" (0, Bz*x, 0)") << endl;
    ss << mPrefix << " eps: " << endl << meps;
    ss << mPrefix << " t01x: " << endl << mt01x;
    ss << mPrefix << " t01y: " << endl << mt01y;

    return ss.str(); 
};

void GrapheneKpParams::update(){
    if(!(is_finite(mdtol) && is_finite(mgamma) && is_finite(ma) 
            && is_finite(mK))){
        throw runtime_error("GrapheneKpParams: invalid k.p parameters.");    
    }
    double Kx = mK, Ky = mK, ax = ma, ay = ma;
    meps = -mgamma*(Kx/ax + Ky/ay)*sz();
    mt01x = (mgamma/(2*ax))*sx()*i + (Kx*mgamma/(2*ax))*sz();
    mt10x = trans(mt01x);
    mt01y = (mgamma/(2*ay))*sy()*i + (Ky*mgamma/(2*ay))*sz();
    mt10y = trans(mt01y);                
}

cxmat GrapheneKpParams::twoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj) const
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1) {
        throw invalid_argument("GrapheneHamGen(): atomi and atomj must should contain one atom each.");
    }
    
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat hmat =  zeros<cxmat>(noi, noj);

    // calculate distance between atom i and atom j
    double xi = atomi.X(0), yi = atomi.Y(0), zi = atomi.Z(0);
    double xj = atomj.X(0), yj = atomj.Y(0), zj = atomj.Z(0);
    
    double dx = abs(xi - xj), dy = abs(yi - yj), dz = abs(zi - zj);           
    double d = sqrt(dx*dx + dy*dy + dz*dz);

    // Assign the the matrix elements based on the distance between 
    // the lattice points.
    if (atomi.Symbol(0) == "D" && atomj.Symbol(0) == "D"){
        // site energy
        if (d <= mdtol){
            hmat = meps;
        // nearest neighbor in x
        }else if(abs(d - ma) <= mdtol && abs(dx - ma) <= mdtol){ 
            // add magnetic field
            dcmplx phase = calcPeierlsPhase(xi, yi, xj, yj);
            
            if (xi > xj){
                hmat = mt10x*phase;
            } else {
                hmat = mt01x*phase;
            }
        //nearest neighbor y
        }else if (abs(d - ma) <= mdtol && abs(dy - ma) <= mdtol){
            // add magnetic field
            dcmplx phase = calcPeierlsPhase(xi, yi, xj, yj);
            
            if(yi > yj){
                hmat = mt10y*phase;
            }else{
                hmat = mt01y*phase;
            }
        }            
    }  
    return hmat;
};

cxmat GrapheneKpParams::twoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj) const
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1) {
        throw invalid_argument("GrapheneHamGen(): atomi and atomj must should contain one atom each.");
    }
    
    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat smat =  zeros<cxmat>(noi, noj);

    // calculate distance between atom i and atom j
    double xi = atomi.X(0), yi = atomi.Y(0), zi = atomi.Z(0);
    double xj = atomj.X(0), yj = atomj.Y(0), zj = atomj.Z(0);
    
    double dx = abs(xi - xj), dy = abs(yi - yj), dz = abs(zi - zj);           
    double d = sqrt(dx*dx + dy*dy + dz*dz);

    // Assign the the matrix elements based on the distance between 
    // the lattice points.
    if (atomi.Symbol(0) == "D" && atomj.Symbol(0) == "D"){
        // overlap matrix of atom i.
        if (d <= mdtol){
            smat = mI;
        }
    }
    return smat;
};


}
}


