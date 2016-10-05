/* 
 * File:   tikp4.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 28, 2014, 08:31 PM
 */

#include "hamiltonian/kp/tikp4.h"

namespace quest{
namespace hamiltonian{


TISurfKpParams4::TISurfKpParams4(const string &prefix):cxhamparams(prefix)
{
    mTitle = "Topological Insulator surface k.p parameters with four spins";   
    mI = eye<cxmat>(4,4);
    
    setDefaultParams();
    update();
}

string TISurfKpParams4::toString() const { 
    stringstream ss;
    ss << cxhamparams::toString() << ":" << endl;
    ss << mPrefix << " a    = " << ma << endl;
    ss << mPrefix << " K    = " << mK << endl;
    ss << mPrefix << " C    = " << mC << endl;
    ss << mPrefix << " A2   = " << mA2 << endl;
    ss << mPrefix << " dtol = " << mdtol << endl;    
    ss << mPrefix << " eps: " << endl << meps;
    ss << mPrefix << " t01x: " << endl << mt01x;
    ss << mPrefix << " t01y: " << endl << mt01y;

    return ss.str(); 
}   

void TISurfKpParams4::update(){

    if(!(is_finite(mdtol) && is_finite(mC) && is_finite(mA2) && is_finite(ma)
            && is_finite(mK))){
        throw runtime_error("TISufrParams: invalid TI parameters.");            
    }
    // To avoid the famous Fermion doubling and its side effects, we are
    // adding sigma_z(kx^2 + ky^2) to the upper 2x2 and subtract from the 
    // lower 2x2. 
    double Kx = mK, Ky = mK, ax = ma, ay = ma;
    cxmat I22 = eye<cxmat>(2,2); 
    cxmat tmp1, tmp2;

    meps = zeros<cxmat>(4,4);
    tmp1 = mC*I22 - mA2*(Kx/ax + Ky/ay)*sz();
    tmp2 = mC*I22 + mA2*(Kx/ax + Ky/ay)*sz();
    meps(span(0,1), span(0,1)) = tmp1;
    meps(span(2,3), span(2,3)) = tmp2;

    mt01x = zeros<cxmat>(4,4);
    tmp1 = (mA2/(2*ax))*sy()*i + (Kx*mA2/(2*ax))*sz();
    tmp2 = (mA2/(2*ax))*sy()*i - (Kx*mA2/(2*ax))*sz();
    mt01x(span(0,1), span(0,1)) = tmp1;
    mt01x(span(2,3), span(2,3)) = tmp2;
    mt10x = trans(mt01x);

    mt01y = zeros<cxmat>(4,4);
    tmp1 = (-mA2/(2*ay))*sx()*i + (Ky*mA2/(2*ay))*sz();
    tmp1 = (-mA2/(2*ay))*sx()*i - (Ky*mA2/(2*ay))*sz();
    mt01y(span(0,1), span(0,1)) = tmp1;
    mt01y(span(2,3), span(2,3)) = tmp2;
    mt10y = trans(mt01y);
}  

cxmat TISurfKpParams4::twoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj) const
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
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
            if (xi > xj){
                hmat = mt10x;
            }else{
                hmat = mt01x;
            }
        //nearest neighbor y
        }else if (abs(d - ma) <= mdtol && abs(dy - ma) <= mdtol){
            if(yi > yj){
                hmat = mt10y;
            }else{
                hmat = mt01y;
            }
        }            
    }
    
    return hmat;
};

cxmat TISurfKpParams4::twoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj) const
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
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
        // site energy
        if (d <= mdtol){
            smat = mI;
        }
    }   
    return smat;
};

}
}




