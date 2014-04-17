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

struct TISurfKpParams: public HamParams{
    //!< k.p parameters
    double ax;           // discretization length in x direction
    double ay;           // discretization length in y direction
    double C;            // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double A2;           // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double Rx;           // parameter to avoid fermion doubling
    double Ry;           // parameter to avoid fermion doubling
    
    //!< tight binding parameters
    cxmat22 I;
    cxmat22 eps;
    cxmat22 t01x;
    cxmat22 t10x;
    cxmat22 t01y;
    cxmat22 t10y;
    
    TISurfKpParams(const string &prefix = ""):HamParams(prefix){
        mTitle = "Topological Insulator surface k.p parameters";
        // default parameters
        I = eye<cxmat>(2,2);    
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
        eps = C*I - 2*A2*(Rx/(ax*ax) + Ry/(ay*ay))*sz();
        t01x =  (A2/(2*ax))*sy()*i + (Rx*A2/(ax*ax))*sz();
        t10x = trans(t01x);
        t01y = (-A2/(2*ay))*sx()*i + (Ry*A2/(ay*ay))*sz();
        t10y = trans(t01y);
    };

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

