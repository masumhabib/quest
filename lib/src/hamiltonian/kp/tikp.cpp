/* 
 * File:   ti.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#include "hamiltonian/kp/tikp.h"

namespace quest{
namespace hamiltonian{


TISurfKpParams::TISurfKpParams(const string &prefix):cxhamparams(prefix)
{
    mTitle = "Topological Insulator surface k.p parameters";
    mI = eye<cxmat>(2,2);
   
    setDefaultParams();
    update();
}
        
string TISurfKpParams::toString() const { 
    stringstream ss;
    ss << cxhamparams::toString() << ":" << endl;
    ss << mPrefix << " a    = " << ma << endl;
    ss << mPrefix << " K    = " << mK << endl;
    ss << mPrefix << " C    = " << mC << endl;
    ss << mPrefix << " A2   = " << mA2 << endl;
    ss << mPrefix << " dtol = " << mdtol << endl;
    ss << mPrefix << " Bz   = " << mBz << endl;
    ss << mPrefix << " A    = " << (mBzGauge == coord::X?" (-Bz*y, 0, 0)":" (0, Bz*x, 0)") << endl;    
    ss << mPrefix << " eps: " << endl << meps;
    ss << mPrefix << " t01x: " << endl << mt01x;
    ss << mPrefix << " t01y: " << endl << mt01y;

    return ss.str(); 
};
    
void TISurfKpParams::update(){

    if(!(is_finite(mdtol) && is_finite(mC) && is_finite(mA2) && is_finite(ma)
             && is_finite(mK))){
        throw runtime_error("TISurfKpParams: invalid TI parameters.");            
    }
    double Kx = mK, Ky = mK, ax = ma, ay = ma;
    meps = mC*mI - mA2*(Kx/ax + Ky/ay)*sz();
    mt01x =  (mA2/(2*ax))*sy()*i + (Kx*mA2/(2*ax))*sz();
    mt10x = trans(mt01x);
    mt01y = (-mA2/(2*ay))*sx()*i + (Ky*mA2/(2*ay))*sz();
    mt10y = trans(mt01y);
} 

cxmat TISurfKpParams::twoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)  const
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
            // add magnetic field
            dcmplx phase = dcmplx(1,0);
            if(abs(mBz) > mBzTol){
                if (mBzGauge == coord::X){ // for A = (-Bz*y, 0, 0)
                    double phi = mfactor*mBz*(xi - xj)*(yi+yj);
                    phase = exp(i*phi);
                }else if (mBzGauge == coord::Y){ // for A = (0, Bz*x, 0)
                    double phi = mfactor*mBz*(yj - yi)*(xi+xj);
                    phase = exp(i*phi);
                }
            }
            
            if (xi > xj){
                hmat = mt10x*phase;
            }else{
                hmat = mt01x*phase;
            }
        //nearest neighbor y
        }else if (abs(d - ma) <= mdtol && abs(dy - ma) <= mdtol){
            // add magnetic field
            dcmplx phase = dcmplx(1,0);
            if(abs(mBz) > mBzTol){
                if (mBzGauge == coord::Y){ // for A = (0, Bz*x, 0)
                    double phi = mfactor*mBz*(yj - yi)*(xi+xj);
                    phase = exp(i*phi);
                }else if (mBzGauge == coord::X){ // for A = (-Bz*y, 0, 0)
                    double phi = mfactor*mBz*(xi - xj)*(yi+yj);
                    phase = exp(i*phi);
                }
            }
            
            if(yi > yj){
                hmat = mt10y*phase;
            }else{
                hmat = mt01y*phase;
            }
        }            
    }
    
    return hmat;
};

cxmat TISurfKpParams::twoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj)  const
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




