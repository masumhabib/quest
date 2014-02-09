/* 
 * File:   kpParams.hpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 4, 2014, 2:38 PM
 */

#ifndef KPPARAMS_HPP
#define	KPPARAMS_HPP


#include <armadillo>

#include "../../maths/constants.h"
#include "../../atoms/Atoms.h"
#include "../../utils/Printable.hpp"

#include "../hamParams.hpp"
#include "../genHam.hpp"

namespace qmicad{
using namespace maths::armadillo;
using namespace maths::constants;


struct KpParams: public HamParams{
    // k.p parameters
    double ax;           // discretization length in x direction
    double ay;           // discretization length in y direction
    
    // tight binding parameters
    cxmat22 I;
    cxmat22 eps;
    cxmat22 t01x;
    cxmat22 t10x;
    cxmat22 t01y;
    cxmat22 t10y;
    
    KpParams(const string &prefix = ""):HamParams(prefix){
        mTitle = "k.p paramters";
        // default parameters          
    }
    
    // Updates internal tight binding parameters calculated using 
    // k.p model. Call it after changing any of the k.p parameters.
    virtual void update(){}; 
    virtual string toString() const {};
        
protected:

};
}

#endif	/* KPPARAMS_HPP */

