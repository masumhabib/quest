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

namespace qmicad{
namespace python{

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


};

}
}
#endif	/* PYLINEARPOT_H */

