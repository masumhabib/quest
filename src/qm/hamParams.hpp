/* 
 * File:   hamParams.hpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 5, 2014, 9:50 AM
 * 
 * Base class for Hamiltonian parameters
 * 
 */

#ifndef HAMPARAMS_HPP
#define	HAMPARAMS_HPP

#include "../atoms/Atoms.h"
#include "../utils/Qbase.hpp"

namespace qmicad{
using namespace constants;


struct HamParams: public Qbase{
    // Parameters required for all Hamiltonin
    double dtol;         // distance tolerance    
    // our periodic table
    vector<Atom> PeriodicTable;
    
    HamParams(const string &prefix = ""):Qbase(" " + prefix){
        mTitle = "Hamiltonian paramters";
        // default parameters          
    }
    
    // Updates internal parameters. Call it after changing any of the 
    // public parameters.
    virtual void update(){}; 
    virtual string toString() const {};
        
protected:

};

}
#endif	/* HAMPARAMS_HPP */

