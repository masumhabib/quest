/* 
 * File:   graphenekp.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#ifndef GRAPHENEKP_H
#define	GRAPHENEKP_H

#include <armadillo>

#include "../../maths/constants.h"
#include "../../atoms/AtomicStruct.h"
#include "../../utils/Printable.hpp"
#include "../hamiltonian.hpp"

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

struct GrapheneKpParams: public HamParams{
    //!< k.p parameters
    double ax;           // discretization length in x direction
    double ay;           // discretization length in y direction
    double gamma;        // h_bar*v_F for graphene
    double Rx;           // parameter to avoid fermion doubling
    double Ry;           // parameter to avoid fermion doubling

    //!< tight binding parameters
    cxmat I;
    cxmat eps;
    cxmat t01x;
    cxmat t10x;
    cxmat t01y;
    cxmat t10y;
    
    GrapheneKpParams(const string &prefix = ""):HamParams(prefix){
        mTitle = "Graphene k.p parameters";
        // default parameters          
        I = eye<cxmat>(2,2);    
    }
    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update(){
        if(!(is_finite(dtol) && is_finite(gamma) && is_finite(ax) 
                && is_finite(ay) && is_finite(Rx) && is_finite(Ry))){
            throw runtime_error("TISufrParams: invalid TI parameters.");    
        }

        computeTBParams();
    }
    
    virtual string toString() const { 
        stringstream ss;
        ss << HamParams::toString() << ":" << endl;
        ss << mPrefix << " dtol  = " << dtol << endl;
        ss << mPrefix << " ax    = " << ax << endl;
        ss << mPrefix << " ay    = " << ay << endl;
        ss << mPrefix << " Rx    = " << Rx << endl;
        ss << mPrefix << " Ry    = " << Rx << endl;
        ss << mPrefix << " gamma = " << gamma << endl;
        ss << mPrefix << " eps: " << endl << eps;
        ss << mPrefix << " t01x: " << endl << t01x;
        ss << mPrefix << " t01y: " << endl << t01y;

        return ss.str(); 
    };  
    
protected:
    void computeTBParams(){
        eps = -2*gamma*(Rx/(ax*ax) + Ry/(ay*ay))*sz();
        t01x = (gamma/(2*ax))*sx()*i + (Rx*gamma/(ax*ax))*sz();
        t10x = trans(t01x);
        t01y = (gamma/(2*ay))*sy()*i + (Ry*gamma/(ay*ay))*sz();
        t10y = trans(t01y);
    };

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

