/* 
 * File:   TI3DKpParams.cpp
 * Author: K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on September 12, 2014, 9:49 AM
 */

#include "hamiltonian/kp/TI3DKpParams.h"

namespace qmicad{namespace hamiltonian{
TI3DKpParams::TI3DKpParams(const string material, 
        const string &prefix):cxhamparams(prefix) 
{
    mTitle = "Topological Insulator four band k.p parameters";   
    mI = eye<cxmat>(4,4);
    mU = zeros<cxmat>(4,4);
    mU(0, 0) = 1;
    mU(1, 2) = 1;
    mU(2, 1) = 1;
    mU(3, 3) = 1;
    
    setDefaultParams();
    if (material == "Bi2Se3"){
        setBi2Se3Params();
    }else{
        throw runtime_error(string("In TI3DKpParams(): ") + material + " is not implemented yet." );
    }
    
    update();
    
}

string TI3DKpParams::toString() const { 
    stringstream ss;
    ss << cxhamparams::toString() << ":" << endl;
    ss << mPrefix << " a    = " << ma << endl;    
    ss << mPrefix << " dtol = " << mdtol << endl;    
    ss << mPrefix << " A1   = " << mA1 << endl;
    ss << mPrefix << " A2   = " << mA2 << endl;
    ss << mPrefix << " B1   = " << mB1 << endl;
    ss << mPrefix << " B2   = " << mB2 << endl;
    ss << mPrefix << " C    = " << mC << endl;
    ss << mPrefix << " D1   = " << mD1 << endl;
    ss << mPrefix << " D2   = " << mD2 << endl;
    ss << mPrefix << " M    = " << mC << endl;
    ss << mPrefix << " eps: " << endl << meps;
    ss << mPrefix << " t01x: " << endl << mt01x;
    ss << mPrefix << " t01y: " << endl << mt01y;
    ss << mPrefix << " t01z: " << endl << mt01z;

    return ss.str(); 
}   

void TI3DKpParams::update(){

    if(!(is_finite(mdtol) && is_finite(ma)
            && is_finite(mA1) && is_finite(mA2)
            && is_finite(mB1) && is_finite(mB2)
            && is_finite(mC)
            && is_finite(mD1) && is_finite(mD2)
            && is_finite(mM)
            )){
        throw runtime_error("In TI3DKpParams::update(): invalid TI parameters.");            
    }
    double ax = ma, ay = ma, az = ma;
    
    cxmat Go(4,4,fill::zeros);
    cxmat Gx(4,4,fill::zeros), Gx2(4,4,fill::zeros);
    cxmat Gy(4,4,fill::zeros), Gy2(4,4,fill::zeros);
    cxmat Gz(4,4,fill::zeros), Gz2(4,4,fill::zeros);
    
    // see TI notes, page 11.
    Gx(0,3) = mA2;
    Gx(1,2) = mA2;
    Gx(2,1) = mA2;
    Gx(3,0) = mA2;
    
    Gy(0,3) = -i*mA2;
    Gy(1,2) = -i*mA2;
    Gy(2,1) = i*mA2;
    Gy(3,0) = i*mA2;

    Gx2(0,0) = mD2 - mB2;
    Gx2(1,1) = mD2 + mB2;
    Gx2(2,2) = mD2 - mB2;
    Gx2(3,3) = mD2 + mB2;
    
    Gy2(0,0) = mD2 - mB2;
    Gy2(1,1) = mD2 + mB2;
    Gy2(2,2) = mD2 - mB2;
    Gy2(3,3) = mD2 + mB2;
    
    Gz(0,1) = mA1;
    Gz(1,0) = mA1;
    Gz(2,3) = -mA1;
    Gz(3,2) = -mA1;

    Gz2(0,0) = mD1 - mB1;
    Gz2(1,1) = mD1 + mB1;
    Gz2(2,2) = mD1 - mB1;
    Gz2(3,3) = mD1 + mB1;
    
    Go(0,0) = mC + mM;
    Go(1,1) = mC - mM;
    Go(2,2) = mC + mM;
    Go(3,3) = mC - mM;
   
    meps = Go + (2/(ax*ax))*Gx2 + (2/(ay*ay))*Gy2 + (2/(az*az))*Gz2;
    mt01x = (-i/(2*ax))*Gx + (1/(ax*ax))*Gx2;
    mt01y = (-i/(2*ay))*Gy + (1/(ay*ay))*Gy2;
    mt01z = (-i/(2*az))*Gz + (1/(az*az))*Gz2;
    
    // change basis to |1, up> |1, dn> |2, up> |2, dn>
    mt01x = trans(mU)*mt01x*mU;
    mt01y = trans(mU)*mt01y*mU;
    mt01z = trans(mU)*mt01z*mU;
    
    mt10x = trans(mt01x);
    mt10y = trans(mt01y);
    mt10z = trans(mt01z);

}  

cxmat TI3DKpParams::twoAtomHam(const AtomicStruct& atomi, 
        const AtomicStruct& atomj) const
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
    }

    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat hmat =  zeros<cxmat>(noi, noj);

    // calculate distance between atom i and atom j
    double xi = atomi.X(0);
    double yi = atomi.Y(0);
    double zi = atomi.Z(0);
    
    double xj = atomj.X(0);
    double yj = atomj.Y(0);
    double zj = atomj.Z(0);

    double dx = abs(xi - xj);
    double dy = abs(yi - yj);           
    double dz = abs(zi - zj);           
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
        //nearest neighbor z
        }else if (abs(d - ma) <= mdtol && abs(dz - ma) <= mdtol){
            if(zi > zj){
                hmat = mt10z;
            }else{
                hmat = mt01z;
            }
        }            

    }
    
    return hmat;
}

cxmat TI3DKpParams::twoAtomOvl(const AtomicStruct& atomi, 
        const AtomicStruct& atomj) const
{
    if (atomi.NumOfAtoms() > 1 || atomj.NumOfAtoms() > 1){
        throw invalid_argument("TISurfHamGen(): atomi and atomj must should contain one atom each.");
    }

    int noi = atomi.NumOfOrbitals();
    int noj = atomj.NumOfOrbitals();

    cxmat smat =  zeros<cxmat>(noi, noj);

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
        if (d <= mdtol){
            smat = mI;
        }
    }   
    return smat;
};



}}

