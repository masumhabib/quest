/* 
 * File:   ti.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#ifndef TIKP_H
#define	TIKP_H

#include <armadillo>

#include "../../utils/mymath.h"
#include "../../utils/Printable.hpp"
#include "../../atoms/Atoms.h"
#include "../genHam.hpp"
#include "kpParams.hpp"

using arma::cx_mat22;
using arma::cx_mat;
using namespace constants;

struct TISurfKpParams: public KpParams{
    // k.p parameters
    double C;            // C parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double A2;           // A2 parameter, see Rev. Mod. Phys. 83, 1057 (2011)
    double K;            // parameter to avoid fermion doubling
        
    TISurfKpParams(const string &prefix = ""):KpParams(prefix){
        mTitle = "Topological Insulator Surface";
        // default parameters
          
        I = eye<cx_mat>(2,2);    
        // Discretized lattice points for k.p Hamiltonian
        PeriodicTable.push_back(Atom(0,  "D",  2, 2));

    }
    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update(){
        
        if(!(is_finite(dtol) && is_finite(C) && is_finite(A2) && is_finite(ax)
                && is_finite(ay) && is_finite(K))){
            throw runtime_error("TISufrParams: invalid TI parameters.");
            
        }
        
        computeTBParams();
    }
    
    virtual string toString() const { 
        stringstream ss;
        ss << mPrefix << " dtol = " << dtol << endl;
        ss << mPrefix << " ax = " << ax << ", ay = " << ay << endl;
        ss << mPrefix << " K = " << K << endl;
        ss << mPrefix << " C = " << C << " A2 = " << A2 << endl;
        ss << mPrefix << " eps = " << endl << eps << endl;
        ss << mPrefix << " t01x " << endl << t01x << endl;
        ss << mPrefix << " t10x " << endl << t10x << endl;
        ss << mPrefix << " t01y " << endl << t01y << endl;
        ss << mPrefix << " t10y " << endl << t10y << endl;

        return ss.str(); 
    };
        
protected:
    void computeTBParams(){
        eps = C*I + (K*A2/ax + K*A2/ay)*sz();
        t01x =  (A2/(2*ax))*sy()*i - (K*A2/(2*ax))*sz();
        t10x = trans(t01x);
        t01y = (-A2/(2*ay))*sx()*i - (K*A2/(2*ay))*sz();
        t10y = trans(t01y);
    };

};
/* Tight binding logic for graphene*/
class TISurfHamGen: public HamGenerator<cx_mat>{
// Fields
protected:    
    TISurfKpParams p;
// Methods
public:
    TISurfHamGen(const TISurfKpParams& p){this->p = p;};
    virtual ~TISurfHamGen(){};
protected:
    virtual cx_mat operator()(const Atoms& atomi, const Atoms& atomj) const;
    
private:
    
};

#endif	/* TIKP_H */

