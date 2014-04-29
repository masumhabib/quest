/* 
 * File:   tikp4.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 28, 2014, 08:31 PM
 */

#ifndef TIKP4_H
#define	TIKP4_H

#include <armadillo>

#include "../../maths/constants.h"
#include "../../utils/Printable.hpp"
#include "../../atoms/AtomicStruct.h"
#include "../hamiltonian.hpp"

namespace qmicad{
namespace hamiltonian{

using boost::static_pointer_cast;
using namespace maths::armadillo;
using namespace maths::constants;

/**
 * TI surface k.p Hamiltonian parameters in four spins basis set.
 * In four spin basis [z_up, z_dn, z_up, z_dn], the hamiltonian 
 * for one site has the form
 *  H_uu  Hud  0    0
 *  H_du  Hdd  0    0
 *  0     0    Huu  Hud
 *  0     0    Hdu  Hdd
 * A quantity, sigma_z*(kx^2 + ky^2) is added to the upper 2x2 and subtracted 
 * from the lowe 2x2 to avoid fermion doubling and to avoid all the side
 * effects that can be introduced by the term sigma_z*(kx^2 + ky^2).
 */
struct TISurfKpParams4: public HamParams{
    //!< k.p parameters
    double ax;           // discretization length in x direction
    double ay;           // discretization length in y direction
    double C;            // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double A2;           // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double Rx;           // parameter to avoid fermion doubling
    double Ry;           // parameter to avoid fermion doubling
    
    //!< tight binding parameters
    cxmat I;
    cxmat eps;
    cxmat t01x;
    cxmat t10x;
    cxmat t01y;
    cxmat t10y;
    
    TISurfKpParams4(const string &prefix = ""):HamParams(prefix){
        mTitle = "Topological Insulator surface k.p parameters with four spins";   
        I = eye<cxmat>(4,4); 
    }
    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update(){
        
        if(!(is_finite(dtol) && is_finite(C) && is_finite(A2) && is_finite(ax)
                && is_finite(ay) && is_finite(Rx) && is_finite(Ry))){
            throw runtime_error("TISufrParams: invalid TI parameters.");
            
        }
        
        computeTBParams();
    }
    
    virtual string toString() const { 
        stringstream ss;
        ss << HamParams::toString() << ":" << endl;
        ss << mPrefix << " dtol = " << dtol << endl;
        ss << mPrefix << " ax   = " << ax << endl;
        ss << mPrefix << " ay   = " << ay << endl;
        ss << mPrefix << " Rx   = " << Rx << endl;
        ss << mPrefix << " Ry   = " << Ry << endl;
        ss << mPrefix << " C    = " << C << endl;
        ss << mPrefix << " A2   = " << A2 << endl;
        ss << mPrefix << " eps: " << endl << eps;
        ss << mPrefix << " t01x: " << endl << t01x;
        ss << mPrefix << " t01y: " << endl << t01y;

        return ss.str(); 
    };
        
protected:
    void computeTBParams(){
        // To avoid the famous Fermion doubling and its side effects, we are
        // adding sigma_z(kx^2 + ky^2) to the upper 2x2 and subtract from the 
        // lower 2x2. 
        cxmat I22 = eye<cxmat>(2,2); 
        cxmat tmp1, tmp2;
        
        eps = zeros<cxmat>(4,4);
        tmp1 = C*I22 - 2*A2*(Rx/(ax*ax) + Ry/(ay*ay))*sz();
        tmp2 = C*I22 + 2*A2*(Rx/(ax*ax) + Ry/(ay*ay))*sz();
        eps(span(0,1), span(0,1)) = tmp1;
        eps(span(2,3), span(2,3)) = tmp2;
        
        t01x = zeros<cxmat>(4,4);
        tmp1 = (A2/(2*ax))*sy()*i + (Rx*A2/(ax*ax))*sz();
        tmp2 = (A2/(2*ax))*sy()*i - (Rx*A2/(ax*ax))*sz();
        t01x(span(0,1), span(0,1)) = tmp1;
        t01x(span(2,3), span(2,3)) = tmp2;
        t10x = trans(t01x);
        
        t01y = zeros<cxmat>(4,4);
        tmp1 = (-A2/(2*ay))*sx()*i + (Ry*A2/(ay*ay))*sz();
        tmp1 = (-A2/(2*ay))*sx()*i - (Ry*A2/(ay*ay))*sz();
        t01y(span(0,1), span(0,1)) = tmp1;
        t01y(span(2,3), span(2,3)) = tmp2;
        t10y = trans(t01y);
    };

};

/** 
 * Tight binding Hamiltonian and overlap matrix for TI surface using 
 *  method k.p.
 */
class TISurfKpHam4: public cxham{
protected:    
public:
    TISurfKpHam4(const TISurfKpParams4& p){
        mhp = shared_ptr<TISurfKpParams4> (new TISurfKpParams4(p));
    };
    virtual ~TISurfKpHam4(){};
protected:
    //!< Generate Hamiltonian between two atoms.
    virtual cxmat genTwoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj);    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat genTwoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj);    
        
};

}
}
#endif	/* TIKP4_H */

