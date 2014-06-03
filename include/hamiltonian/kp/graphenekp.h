/* 
 * File:   graphenekp.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#ifndef GRAPHENEKP_H
#define	GRAPHENEKP_H

#include <armadillo>

#include "maths/constants.h"
#include "atoms/AtomicStruct.h"
#include "utils/Printable.hpp"
#include "hamiltonian/hamiltonian.hpp"

namespace qmicad{
namespace hamiltonian{

using boost::static_pointer_cast;
using namespace maths::armadillo;
using namespace maths::constants;

class GrapheneKpHam;
/**
 * Graphene k.p Hamiltonian parameters with two (pseudo-)spins per site.
 * In two spin basis [z_up, z_dn], the hamiltonian 
 * for one site has the form
 *  H_uu  Hud
 *  H_du  Hdd
 * A quantity, sigma_z*(kx^2 + ky^2) is added to the hamiltonian to avoid 
 * the fermion doubling problem. 
 * See: Phys. Rev. B 86, 085131 (2012) and the references therein.
 */
class GrapheneKpParams: public HamParams{
    friend class GrapheneKpHam;
protected:
    //!< k.p parameters
    double ma;            // discretization length in x direction
    double mgamma;        // h_bar*v_F for graphene
    double mK;            // parameter to avoid fermion doubling

    //!< tight binding parameters
    cxmat mI;
    cxmat meps;
    cxmat mt01x;
    cxmat mt10x;
    cxmat mt01y;
    cxmat mt10y;
    
public:
    GrapheneKpParams(const string &prefix = ""):HamParams(prefix){
        mTitle = "Graphene k.p parameters";
        // default parameters          
        mI = eye<cxmat>(2,2);    
    }
        
    double a(){return ma; }
    void   a(double newa ){ ma = newa; update(); }
    double gamma(){return mgamma; }
    void   gamma(double newgamma ){ mgamma = newgamma; update(); }
    double K(){return mK; }
    void   K(double newK ){ mK = newK; update(); }

    
    virtual string toString() const { 
        stringstream ss;
        ss << HamParams::toString() << ":" << endl;
        ss << mPrefix << " dtol  = " << mdtol << endl;
        ss << mPrefix << " a     = " << ma << endl;
        ss << mPrefix << " K    = " << mK  << endl;
        ss << mPrefix << " gamma = " << mgamma << endl;
        ss << mPrefix << " eps: " << endl << meps;
        ss << mPrefix << " t01x: " << endl << mt01x;
        ss << mPrefix << " t01y: " << endl << mt01y;

        return ss.str(); 
    };
    
protected:
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update(){
        if(!(is_finite(mdtol) && is_finite(mgamma) && is_finite(ma) 
                && is_finite(mK))){
            throw runtime_error("TISufrParams: invalid TI parameters.");    
        }
        double Kx = mK, Ky = mK, ax = ma, ay = ma;
        meps = -mgamma*(Kx/ax + Ky/ay)*sz();
        mt01x = (mgamma/(2*ax))*sx()*i + (Kx*mgamma/(2*ax))*sz();
        mt10x = trans(mt01x);
        mt01y = (mgamma/(2*ay))*sy()*i + (Ky*mgamma/(2*ay))*sz();
        mt10y = trans(mt01y);                
    }    
};

/** 
 * Tight binding Hamiltonian and overlap matrix for graphene using 
 *  method k.p.
 */
class GrapheneKpHam: public cxham{
protected:    
public:
    GrapheneKpHam(const GrapheneKpParams& p){
        mhp = shared_ptr<GrapheneKpParams> (new GrapheneKpParams(p));
    };
    virtual ~GrapheneKpHam(){};
    
protected:
    //!< Generate Hamiltonian between two atoms.
    virtual cxmat genTwoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj);    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat genTwoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj);

};

}
}
#endif	/* GRAPHENEKP_H */

