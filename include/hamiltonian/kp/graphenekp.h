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
class GrapheneKpParams: public cxhamparams{
public:
    GrapheneKpParams(const string &prefix = "");
    
    double a() const {return ma; }
    void   a(double newa ){ ma = newa; update(); }
    double gamma() const {return mgamma; }
    void   gamma(double newgamma ){ mgamma = newgamma; update(); }
    double K() const {return mK; }
    void   K(double newK ){ mK = newK; update(); }

    
    virtual string toString() const;

    //!< Generate Hamiltonian between two atoms.
    virtual cxmat twoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj) const;    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat twoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj) const;
    
private:
    //!< Default parameters.
    void setDefaultParams(){
        ma       = 10.0, 
        mK       = 1.165, 
        mgamma   = 3.16*1.42*3/2, 
        mdtol    = 1E-3; 
        mortho   = true;
        mpt.add(0, "D", 2, 2);
        mBz      = 0;
        mfactor  = 1E-20*q/hbar/2;
    };        
    //!< Updates internal tight binding parameters calculated using 
    //!< k.p model. Call it after changing any of the k.p parameters.
    virtual void update();

private:
    //!< k.p parameters
    double ma;            //!< discretization length in x direction
    double mK;            //!< parameter to avoid fermion doubling
    double mgamma;        //!< h_bar*v_F for graphene

    //!< calculated tight binding matrices
    cxmat mI;
    cxmat meps;
    cxmat mt01x;
    cxmat mt10x;
    cxmat mt01y;
    cxmat mt10y;
    
    double mfactor;      //!< pre-factor for magnetic phase.
    
};

}
}
#endif	/* GRAPHENEKP_H */

