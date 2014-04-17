/* 
 * File:   constants.h
 * Copyright (C) 2013  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 7, 2013, 10:09 PM
 * 
 * Description: Type definitions and math constants
 * 
 */

#ifndef CONSTANTS_H
#define	CONSTANTS_H

#include "arma.hpp"

namespace maths{
namespace constants{
    
using namespace armadillo;

// cmat22: 2x2 matrix helper class that has a 
// constructor which can be used to initialie a const object;
class cmat22{
public:
    cmat22(dcmplx m00, dcmplx m01, dcmplx m10, dcmplx m11){
        m << m00 << m01 << endr << m10 << m11 << endr;
    }
    
    const cxmat& operator()() const{
        return m;
    }
private:
    cxmat m;
};

    // Math 
    const dcmplx        i       = dcmplx(0,1.0);        // imaginary unit number
    const double        eln     = 2.7182818284;         // base of natural logarithm
    const double        pi      = 3.141592653589793;    // The number pi.
    
    // Physics
    const double        h       = 6.62606957E-34;       // Planck constant in: J.s
    const double        hbar    = 1.054571726E-34;      // Reduced Planck constant: J.s
    const double        c       = 299792458;            // Speed of light: m/s
    const double        mub     = 9.27400968E-24;       // Bohr magneton
    
    const double        e       = 1.60217657E-19;       // charge of an electron in: C
    const double        q       = 1.60217657E-19;       // charge of an electron in: C
    const double        me      = 9.10938291E-31;       // Mass of an electron: kg  

    // Pauli matrices
    const cmat22             sx = cmat22(0, 1, 1, 0);
    const cmat22             sy = cmat22(0,-i, i, 0);
    const cmat22             sz = cmat22(1, 0, 0,-1);
    
};

}

#endif	/* MYMATH_H */

