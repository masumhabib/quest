/* 
 * File:   PyBandStruct.h
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on August 17, 2014, 7:40 PM
 */

#ifndef PYBANDSTRUCT_H
#define	PYBANDSTRUCT_H

#include "utils/std.hpp"
#include "band/BandStruct.h"
#include "boostpython.hpp"
#include "npyarma/npyarma.h"

namespace qmicad{namespace python{

namespace bp = boost::python;
using qmicad::python::npy2mat;
using qmicad::python::npy2col;
using boost::make_shared;
using namespace band;

class PyBandStruct: public BandStruct {
// Methods    
public:
    PyBandStruct(const Workers &workers, uint nn,  bool orthoBasis = true, 
            bool calcEigV = false, const string &prefix = "");
    
    //void    H(bp::object H, int ineigh);
    //void    S(bp::object S, int ineigh); 
    void H(const cxmat& H, int ineigh);
    void S(const cxmat& H, int ineigh);
};

}}
#endif	/* PYBANDSTRUCT_H */

