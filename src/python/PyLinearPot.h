/* 
 * File:   PyLinearPot.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 10, 2014, 10:46 PM
 */

#ifndef PYLINEARPOT_H
#define	PYLINEARPOT_H

#include "../potential/linearPot.h"
#include "PyAtomicStruct.h"

#include <boost/shared_ptr.hpp>

namespace qmicad{
namespace python{
using boost::shared_ptr;

class PyLinearPot: public LinearPot {
public:
    PyLinearPot(const PyAtomicStruct &atoms, const string &prefix = ""):
        LinearPot(atoms, prefix)
    {   
    }
        
    friend ostream& operator << (ostream & out, const PyLinearPot &p){
        out << p.toString();
        return out;
    }

    void VLR(int ilr, double Vl, double Vr, double Vt = 0, double Vb = 0){
        LinearPot::VLR(ilr, Vl, Vr, Vt, Vb);
    }
    
    shared_ptr<vec> toOrbPot(uint start, uint end){
        return LinearPot::toOrbPot(span(start, end));
    }

};

}
}
#endif	/* PYLINEARPOT_H */

