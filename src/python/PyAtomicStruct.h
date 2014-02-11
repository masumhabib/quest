/* 
 * File:   PyAtomicStruct.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 9, 2014, 10:36 AM
 */

#ifndef PYATOMICSTRUCT_H
#define	PYATOMICSTRUCT_H

#include "../atoms/AtomicStruct.h"

namespace qmicad{
namespace python{
namespace armadillo = maths::armadillo;

struct PyPeriodicTable:PeriodicTable{
    void add(uint ia, string sym, uint ne, uint no){
        PeriodicTable::add(ia, sym, ne, no);
    }
};

class PyAtomicStruct: public AtomicStruct{
public:
    PyAtomicStruct(const string& gjfFileName ):AtomicStruct(gjfFileName ){
    };
    
    PyAtomicStruct(const string& gjfFileName, const PyPeriodicTable &pTable):
        AtomicStruct(gjfFileName, pTable)
    {
    }
    
    PyAtomicStruct(uint nl, uint nw, double ax, double ay,
        const PyPeriodicTable &pTable): AtomicStruct(nl, nw, ax, ay, pTable)
    {
    }
    
    PyAtomicStruct span(uint start, uint end){
        AtomicStruct tmp = AtomicStruct::operator ()(armadillo::span(start, end));
        return *(static_cast<PyAtomicStruct*>(&tmp));
    }
    
    friend ostream& operator << (ostream & out, const PyAtomicStruct &p){
        out << p.toString();
        return out;
    }

};

}
}
#endif	/* PYATOMICSTRUCT_H */

