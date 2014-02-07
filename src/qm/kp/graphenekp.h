/* 
 * File:   graphenekp.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#ifndef GRAPHENEKP_H
#define	GRAPHENEKP_H

#include <armadillo>

#include "../../utils/mymath.h"
#include "../../atoms/Atoms.h"
#include "../genHam.hpp"
#include "../../utils/Printable.hpp"

#include "kpParams.hpp"

using arma::cx_mat22;
using arma::cx_mat;
using namespace constants;

struct GrapheneKpParams: public KpParams{
    // k.p parameters
    double gamma;        // h_bar*v_F for graphene
    double K;            // parameter to avoid fermion doubling
        
    GrapheneKpParams(const string &prefix = ""):KpParams(prefix){
        mTitle = "Graphene k.p parameters";
        // default parameters
          
        I = eye<cx_mat>(2,2);    
        // Discretized lattice points for k.p Hamiltonian
        PeriodicTable.push_back(Atom(0,  "D",  2, 2));

    }
    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update(){
        if(!(is_finite(dtol) && is_finite(gamma) && is_finite(ax) 
                && is_finite(ay) && is_finite(K))){
            throw runtime_error("TISufrParams: invalid TI parameters.");    
        }

        computeTBParams();
    }
    
    virtual string toString() const { 
        stringstream ss;
        ss << mPrefix << " dtol = " << dtol << endl;
        ss << mPrefix << " ax = " << ax << endl;
        ss << mPrefix << " ay = " << ay << endl;
        ss << mPrefix << " K = " << K << endl;
        ss << mPrefix << " gamma = " << gamma << endl;
        ss << mPrefix << " eps: " << endl << eps << endl;
        ss << mPrefix << " t01x: " << endl << t01x << endl;
        ss << mPrefix << " t01y: " << endl << t01y << endl;

        return ss.str(); 
    };  
    
protected:
    void computeTBParams(){
        eps = (K*gamma/ax + K*gamma/ay)*sz();
        t01x = (gamma/(2*ax))*sx()*i - (K*gamma/(2*ax))*sz();
        t10x = trans(t01x);
        t01y = (gamma/(2*ay))*sy()*i - (K*gamma/(2*ay))*sz();
        t10y = trans(t01y);
    };

};
/* Tight binding logic for graphene*/
class GrapheneKpHamGen: public HamGenerator<cx_mat>{
// Fields
protected:    
    GrapheneKpParams p;
// Methods
public:
    GrapheneKpHamGen(const GrapheneKpParams& p){this->p = p;};
    virtual ~GrapheneKpHamGen(){};
    virtual cx_mat operator()(const Atoms& atomi, const Atoms& atomj) const;
    
private:

};

#endif	/* GRAPHENEKP_H */

