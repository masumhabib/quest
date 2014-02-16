/* 
 * File:   PyGeometry.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 15, 2014, 12:15 AM
 */

#ifndef PYGEOMETRY_H
#define	PYGEOMETRY_H

#include "../maths/geometry.hpp"
#include <python2.6/Python.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/tuple.hpp>

namespace qmicad{
namespace python{
using namespace maths::geometry;
using boost::python::pickle_suite;
using boost::python::tuple;
using namespace boost::python;

class PyPoint: public point{
public:
    PyPoint(const double &x, const double &y): point(x,y){
    }
    
    PyPoint(const point& p){
        this->set<0>(p.get<0>());
        this->set<1>(p.get<1>());
    }
    
    double x() const {return get<0>(); };
    double y() const {return get<1>(); };
};

class PyQuadrilateral: public SimpleQuadrilateral{
public:
    PyQuadrilateral(const PyPoint &lb, const PyPoint &rb, const PyPoint &rt, 
    const PyPoint &lt): SimpleQuadrilateral(lb, rb, rt, lt){
        
    }
    
    PyPoint plb() const { return lb; };
    PyPoint prb() const { return rb; };
    PyPoint prt() const { return rt; };
    PyPoint plt() const { return lt; };
    
    
};

struct PyPointPickler: public pickle_suite{
    static tuple getinitargs(const PyPoint& p){
        return make_tuple(p.x(), p.y());
    }
}; 

struct PyQuadrilateralPickler: public pickle_suite{
    static tuple getinitargs(const PyQuadrilateral& ql){
        return make_tuple(ql.plb(), ql.prb(), ql.prt(), ql.plt());
    }
}; 

}
}
#endif	/* PYGEOMETRY_H */

