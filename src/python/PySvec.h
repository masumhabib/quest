/* 
 * File:   PySvec.h
 * Author: masum
 *
 * Created on February 17, 2014, 2:04 PM
 */

#ifndef PYSVEC_H
#define	PYSVEC_H

#include "../maths/svec.h"
#include "../utils/vout.h"

namespace qmicad{
namespace python{
using namespace maths::spacevec;
using namespace utils::stds;

struct PySvec: public svec{
public:
    PySvec(double x = 0, double y = 0, double z = 0):svec(){
        setx(x);
        setx(y);
        setx(z);
    }
    
    double getx(){
        dout << dbg << " X = " << X << endl;
        dout << dbg << " x = " << (*this)(X) << endl;
        return (*this)(X);
    }

    double gety(){
        return (*this)(Y);
    }

    double getz(){
        return (*this)(Z);
    }

    void setx(double x){
        (*this)(X) = x;
    }

    void sety(double y){
        (*this)(Y) = y;
    }

    void setz(double z){
        (*this)(Z) = z;
    }
    
};

}
}
#endif	/* PYSVEC_H */

