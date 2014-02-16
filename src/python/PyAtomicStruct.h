/* 
 * File:   PyAtomicStruct.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 9, 2014, 10:36 AM
 */

#ifndef PYATOMICSTRUCT_H
#define	PYATOMICSTRUCT_H

#include "../atoms/AtomicStruct.h"
#include "../utils/std.hpp"

#include <python2.6/Python.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/str.hpp>

namespace qmicad{
namespace python{
using namespace boost::python;
namespace armadillo = maths::armadillo;
using namespace utils::stds;

struct PyAtom: public Atom{
    PyAtom(uint ia, string sym, uint ne, uint no):
        Atom(ia, sym, ne, no)
    {
    }
    
    PyAtom():Atom()
    {
    };
    
    PyAtom(const Atom& a){
        this->ia = a.ia;
        this->sym = a.sym;
        this->ne = a.ne;
        this->no = a.no;
    }
    
};

struct PyPeriodicTable: public PeriodicTable{
    void add(uint ia, string sym, uint ne, uint no){
        PeriodicTable::add(ia, sym, ne, no);
    }
    
    void add(const PyAtom &a){
        PeriodicTable::add(a);
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

struct PyAtomPickler : public pickle_suite{
    static tuple getinitargs(const PyAtom& pa){
        return make_tuple(pa.ia, pa.sym, pa.ne, pa.no);
    }
};

 struct PyPeriodicTablePickler : public pickle_suite{
//    static
//    boost::python::tuple
//    getinitargs(const world& w)
//    {
//        using namespace boost::python;
//        return make_tuple(w.get_country());
//    }

    static object getstate(const PyPeriodicTable& pt){
        PyPeriodicTable::piter it;
        // Convert periodic table to string.
        stringstream out;
        for (it = pt.elements.begin(); it != pt.elements.end(); ++it){
            out << it->second.ia << " " << it->second.sym 
                << " " << it->second.ne << " " << it->second.no << endl;
        }
        
        // Return the string to pickle for storage.
        return str(out.str());
    }

    static void setstate(PyPeriodicTable& pt, object state){
        
        // Extract string from pickled state.
        str s = extract<str> (state)();
        string st = extract<string> (s)();    
        stringstream in(st);
    
        // Create the periodic table.
        PyAtom a;
        while(in >> a.ia && in >> a.sym && in >> a.ne && in >> a.no){
            pt.add(a);
        }           
    }
};

}
}
#endif	/* PYATOMICSTRUCT_H */

