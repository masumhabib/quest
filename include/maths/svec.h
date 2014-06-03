/* 
 * File:   svec.h
 * Copyright (C) 2013  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 7, 2013, 9:51 AM
 * 
 * Description: Spetial vector type definition.
 * 
 */

#ifndef SVECT_H
#define	SVECT_H

#include "maths/arma.hpp"


namespace maths{
namespace spvec{

using armadillo::row3;
using armadillo::row2;

namespace coord{
    // index definition
    static const int X = 0;
    static const int Y = 1;
    static const int Z = 2;  
}
 
typedef row3     svec;  // special vector is a three component row vector
typedef row3     svec3;  // special vector is a three component row vector
typedef row2     svec2; // two dimensional vector
   
/*
 * Python wrapper for svec
 */    
class SVec: public svec{
public:
    SVec(double x = 0, double y = 0, double z = 0):svec(){
        setx(x);
        sety(y);
        setz(z);
    }
    
    
    SVec(const svec & rhs):svec(){
        setx(rhs(coord::X));
        sety(rhs(coord::Y));
        setz(rhs(coord::Z));
    }
    
    // Access functions
    double  getx() const    { return (*this)(coord::X); }
    void    setx(double x)  { (*this)(coord::X) = x; }
    double  gety() const    { return (*this)(coord::Y); }
    void    sety(double y)  { (*this)(coord::Y) = y; }
    double  getz() const    { return (*this)(coord::Z); }
    void    setz(double z)  { (*this)(coord::Z) = z; }
    
    svec2 xy() const{ 
        svec2 r2; 
        r2(coord::X) = getx(); 
        r2(coord::Y) = gety(); 
        return r2; 
    }
    
    void xy(double x, double y)  { setx(x); sety(y);}
    
};

}
}
#endif	/* SVECT_H */

