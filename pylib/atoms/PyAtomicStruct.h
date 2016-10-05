/* 
 * File:   PyAtomicStruct.h
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on August 17, 2014, 10:50 AM
 */

#ifndef PYATOMICSTRUCT_H
#define	PYATOMICSTRUCT_H

#include "boostpython.hpp"
#include "atoms/AtomicStruct.h"

/**
 * Python exporters.
 */
namespace quest{ namespace python{
using namespace atoms;
using namespace utils::stds;
using namespace boost::python;

/**
 * @TODO: Implement pickling either using boost::serialization or using python.
 */

/**
 * For python pickle.
 */

struct AtomPickler : public pickle_suite{
    static tuple getinitargs(const Atom& a){
        return make_tuple(a.ia, a.sym, a.ne, a.no);
    }
};


struct AtomicStructPickler : public pickle_suite{

    static object getstate(const AtomicStruct& as){
        ostringstream os;
        boost::archive::text_oarchive oa(os);
        oa << as;
        // Return the string to pickle for storage.
        return str(os.str());
    }

    static void setstate(AtomicStruct& as, object state){
        str s = extract<str> (state)();
        string st = extract<string> (s)();
        istringstream is (st);

        boost::archive::text_iarchive ia (is);
        ia >> as;
    }
};

 struct PeriodicTablePickler : public pickle_suite{

    static object getstate(const PeriodicTable& pt){
        PeriodicTable::cpiter it;
        // Convert periodic table to string.
        stringstream out;
        for (it = pt.elements.begin(); it != pt.elements.end(); ++it){
            out << it->second.ia << " " << it->second.sym 
                << " " << it->second.ne << " " << it->second.no << endl;
        }
        
        // Return the string to pickle for storage.
        return str(out.str());
    }

    static void setstate(PeriodicTable& pt, object state){
        
        // Extract string from pickled state.
        str s = extract<str> (state)();
        string st = extract<string> (s)();    
        stringstream in(st);
    
        // Create the periodic table.
        Atom a;
        while(in >> a.ia && in >> a.sym && in >> a.ne && in >> a.no){
            pt.add(a);
        }           
    }
};

}}

#endif	/* PYATOMICSTRUCT_H */

