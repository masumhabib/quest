/* 
 * File:   arma.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on March 9, 2014, 4:15 PM
 */

#include "arma.hpp"
#include "../python/boostpython.hpp"


/**
 * Python exporter for our favorite armadillo matrices.
 */
namespace qmicad{
namespace python{


void export_cxmat(){
    using namespace maths::armadillo;
    
    class_<cxmat, shared_ptr<cxmat> >("cxmatp", no_init)
    ;
    
}

void export_mat(){
    using namespace maths::armadillo;
    
    class_<mat, shared_ptr<mat> >("matp", no_init)
    ;

}

void export_vec(){
    using namespace maths::armadillo;
    
    class_<vec, shared_ptr<vec> >("vecp", no_init)
    ;
}


}
}