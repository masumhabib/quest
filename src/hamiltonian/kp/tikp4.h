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

class TISurfKpHam4;

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
class TISurfKpParams4: public HamParams{
    friend class TISurfKpHam4;
protected:
    //!< k.p parameters
    double ma;            // discretization length in y direction
    double mC;            // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double mA2;           // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double mK;            // parameter to avoid fermion doubling
    
    //!< tight binding parameters
    cxmat mI;
    cxmat meps;
    cxmat mt01x;
    cxmat mt10x;
    cxmat mt01y;
    cxmat mt10y;
    
public:
    TISurfKpParams4(const string &prefix = ""):HamParams(prefix){
        mTitle = "Topological Insulator surface k.p parameters with four spins";   
        mI = eye<cxmat>(4,4); 
    }
    double a(){return ma; }
    void   a(double newa ){ ma = newa; update(); }
    double C(){return mC; }
    void   C(double newC ){ mC = newC; update(); }
    double A2(){return mA2; }
    void   A2(double newA2 ){ mA2 = newA2; update(); }
    double K(){return mK; }
    void   K(double newK ){ mK = newK; update(); }
    
    virtual string toString() const { 
        stringstream ss;
        ss << HamParams::toString() << ":" << endl;
        ss << mPrefix << " dtol = " << mdtol << endl;
        ss << mPrefix << " a    = " << ma << endl;
        ss << mPrefix << " K    = " << mK << endl;
        ss << mPrefix << " C    = " << mC << endl;
        ss << mPrefix << " A2   = " << mA2 << endl;
        ss << mPrefix << " eps: " << endl << meps;
        ss << mPrefix << " t01x: " << endl << mt01x;
        ss << mPrefix << " t01y: " << endl << mt01y;

        return ss.str(); 
    };
    
protected:
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update(){
        
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

