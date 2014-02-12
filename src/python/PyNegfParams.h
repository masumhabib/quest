/* 
 * File:   PyNegfParams.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 8, 2014, 2:28 AM
 */

#ifndef PYNEGFPARAMS_H
#define	PYNEGFPARAMS_H

#include "../negf/CohRgfa.h"

#include <boost/shared_ptr.hpp>

namespace qmicad{
namespace python{
using boost::shared_ptr;

struct PyNegfParams:public NegfParams{
public:
    PyNegfParams(uint nb):NegfParams(){
        this->nb = nb;
        // initialize the matrix containers.
        NegfParams::H0.set_size(nb);
        NegfParams::S0.set_size(nb);
        NegfParams::Hl.set_size(nb+1);
        NegfParams::Hl.set_size(nb+1);
        NegfParams::V.set_size(nb);
        
    };
    
    void H0(shared_ptr<cxmat> H0, uint it){
        NegfParams::H0(it) = H0;
    }
    
    void S0(shared_ptr<cxmat> S0, uint it){
        NegfParams::S0(it) = S0;
    }
    
    void Hl(shared_ptr<cxmat> Hl, uint it){
        NegfParams::Hl(it) = Hl;
    }
    
    void Sl(shared_ptr<cxmat> Sl, uint it){
        NegfParams::Sl(it) = Sl;
    }
    
    void V(shared_ptr<vec> V, uint it){
        NegfParams::V(it) = V;
    }
    
    friend ostream& operator << (ostream & out, const PyNegfParams &p){
        out << p.toString();
        return out;
    }
};

}
}
#endif	/* PYNEGFPARAMS_H */

