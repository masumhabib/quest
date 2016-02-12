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
#include "hamiltonian/kp/dirackp.h"
#include "hamiltonian/hamiltonian.hpp"

namespace qmicad{
namespace hamiltonian{

using namespace maths::armadillo;
using namespace maths::constants;

/**
 * Base class for Graphene k.p Hamiltonian parameters.
 */
class GrapheneKpParams: public DiracKpParams {
public:
    GrapheneKpParams(const string &prefix = "") : DiracKpParams(prefix) {
        mTitle  = "Graphene k.p parameters";
    };
    
    double gamma() const {return mgamma; }
    void   gamma(double newgamma ){ mgamma = newgamma; update(); }
    
    virtual string toString() const;
protected:
    //!< Default parameters.
    virtual void setDefaultParams(){
        ma       = 10.0, 
        mK       = 1.165, 
        mgamma   = 3.16*1.42*3/2, 
        mdtol    = 1E-3; 
        mortho   = true;
        mBz      = 0;
    };        
    //!< Updates internal tight binding parameters calculated using 
    //!< k.p model. Call it after changing any of the k.p parameters.
    virtual void update() {};

protected:
    //!< k.p parameters
    double mgamma;        //!< h_bar*v_F for graphene
};

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
class GrapheneOneValleyKpParams: public GrapheneKpParams {
public:
    GrapheneOneValleyKpParams(const string &prefix = "");
    
protected:
    //!< Default parameters.
    virtual void setDefaultParams(){
        GrapheneKpParams::setDefaultParams();
        mpt.add(0, "D", 2, 2);
    }

    //!< Updates internal tight binding parameters calculated using 
    //!< k.p model. Call it after changing any of the k.p parameters.
    virtual void update();

private:
};

/**
 * Graphene k.p Hamiltonian parameters with four (pseudo-)spins per site.
 * In two spin basis [z_up_K, z_dn_K, z_up_K', z_dn_K'], the hamiltonian 
 * for one site has the form
 *  H_uu  Hud  0    0
 *  H_du  Hdd  0    0
 *  0     0    Huu  Hud
 *  0     0    Hdu  Hdd

 * A quantity, sigma_z*(kx^2 + ky^2) is added to the hamiltonian to avoid 
 * the fermion doubling problem. 
 * See: Phys. Rev. B 86, 085131 (2012) and the references therein.
 */
class GrapheneTwoValleyKpParams: public GrapheneKpParams {
public:
    GrapheneTwoValleyKpParams(const string &prefix = "");
    
protected:
    virtual void setDefaultParams(){
        GrapheneKpParams::setDefaultParams();
        mpt.add(0, "D", 4, 4);
    }

    //!< Updates internal tight binding parameters calculated using 
    //!< k.p model. Call it after changing any of the k.p parameters.
    virtual void update();

private:
};

}
}

#endif	/* GRAPHENEKP_H */

