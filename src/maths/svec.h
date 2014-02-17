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

#include "arma.hpp"


namespace maths{
namespace spacevec{
using armadillo::row;
using armadillo::row3;
using armadillo::row2;

    // index definition
    static const int X = 0;
    static const int Y = 1;
    static const int Z = 2;  

    typedef row3     svec;  // special vector is a three component row vector
    typedef row2     xyvec; // two dimensional vector
/*    
class SVec: public row{
public:
    SVec(double x = 0, double y = 0, double z = 0):row(3){
        setx(x);
        setx(y);
        setx(z);
    }
    
    // Access functions
    double  x() const    { return (*this)(X); }
    void    x(double x)  { (*this)(X) = x; }
    double  y() const    { return (*this)(Y); }
    void    y(double y)  { (*this)(Y) = y; }
    double  z() const    { return (*this)(Z); }
    void    z(double z)  { (*this)(Z) = z; }

    //!< For python wrapper.
    double  getx() const  { return x(); };
    void    setx(double x){ x(x);}
    double  gety() const  { return y(); };
    void    sety(double y){ y(y);}
    double  getz() const  { return z(); };
    void    setz(double z){ z(z);}
    
public:
    // index definition
    static const int X = 0;
    static const int Y = 1;
    static const int Z = 2;  
};
   */ 
}
}
#endif	/* SVECT_H */

