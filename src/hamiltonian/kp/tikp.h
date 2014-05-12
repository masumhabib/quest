/* 
 * File:   ti.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#ifndef TIKP_H
#define	TIKP_H

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

class TISurfKpHam;

/**
 * TI surface k.p Hamiltonian parameters with two spins per site.
 * In two spin basis [z_up, z_dn], the hamiltonian 
 * for one site has the form
 *  H_uu  Hud
 *  H_du  Hdd
 * A quantity, sigma_z*(kx^2 + ky^2) is added to the hamiltonian to avoid 
 * the fermion doubling problem. 
 * See: Phys. Rev. B 86, 085131 (2012) and the references therein.
 */
class TISurfKpParams: public HamParams{
friend class TISurfKpHam;
protected:
    //!< k.p parameters
    double ma;           // discretization length in x direction
    double mC;            // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double mA2;           // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double mK;           // parameter to avoid fermion doubling
    
    //!< tight binding parameters
    cxmat mI;
    cxmat meps;
    cxmat mt01x;
    cxmat mt10x;
    cxmat mt01y;
    cxmat mt10y;

public:    
    TISurfKpParams(const string &prefix = ""):HamParams(prefix){
        mTitle = "Topological Insulator surface k.p parameters";
        // default parameters
        mI = eye<cxmat>(2,2);    
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
        double Kx = mK, Ky = mK, ax = ma, ay = ma;
        meps = mC*mI - mA2*(Kx/ax + Ky/ay)*sz();
        mt01x =  (mA2/(2*ax))*sy()*i + (Kx*mA2/(2*ax))*sz();
        mt10x = trans(mt01x);
        mt01y = (-mA2/(2*ay))*sx()*i + (Ky*mA2/(2*ay))*sz();
        mt10y = trans(mt01y);
    }    
};

/** 
 * Tight binding Hamiltonian and overlap matrix for TI surface using 
 *  method k.p.
 */
class TISurfKpHam: public cxham{
protected:    
public:
    TISurfKpHam(const TISurfKpParams& p){
        mhp = shared_ptr<TISurfKpParams> (new TISurfKpParams(p));
    };
    virtual ~TISurfKpHam(){};
protected:
    //!< Generate Hamiltonian between two atoms.
    virtual cxmat genTwoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj);    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat genTwoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj);    
        
};

}
}
#endif	/* TIKP_H */

