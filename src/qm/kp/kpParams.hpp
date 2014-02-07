/* 
 * File:   kpParams.hpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 4, 2014, 2:38 PM
 */

#ifndef KPPARAMS_HPP
#define	KPPARAMS_HPP


#include <armadillo>

#include "../../utils/mymath.h"
#include "../../atoms/Atoms.h"
#include "../../utils/Printable.hpp"

#include "../hamParams.hpp"
#include "../genHam.hpp"

using arma::cx_mat22;
using arma::cx_mat;
using namespace constants;


struct KpParams: public HamParams{
    // k.p parameters
    double ax;           // discretization length in x direction
    double ay;           // discretization length in y direction
    
    // tight binding parameters
    cx_mat22 I;
    cx_mat22 eps;
    cx_mat22 t01x;
    cx_mat22 t10x;
    cx_mat22 t01y;
    cx_mat22 t10y;
    
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


#endif	/* KPPARAMS_HPP */

