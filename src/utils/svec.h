/* 
 * File:   svec.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 7, 2013, 9:51 AM
 * 
 * Description: Spetial vector type definition.
 * 
 */

#ifndef SVECT_H
#define	SVECT_H

#include <armadillo>

namespace spacevec{
    // index definition
    static const int X = 0;
    static const int Y = 1;
    static const int Z = 2;  
    
    typedef arma::rowvec3 svec;  // special vector is a three component row vector
    typedef arma::rowvec2 xyvec; // two dimensional vector
    
}

#endif	/* SVECT_H */

