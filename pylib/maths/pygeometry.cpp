/* 
 * File:   geometry.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 3, 2014, 1:45 PM
 */

#include "maths/geometry.hpp"
#include "boostpython.hpp"


namespace qmicad{
namespace python{

using maths::geometry::point;
using maths::geometry::quadrilateral;

/**
 * Geometry python wrappers.
 */

struct PointPickler: public pickle_suite{
    static tuple getinitargs(const point& p){
        return make_tuple(p.get<0>(), p.get<1>());
    }
}; 

struct QuadrilateralPickler: public pickle_suite{
    static tuple getinitargs(const quadrilateral& ql){
        return make_tuple(ql.lb, ql.rb, ql.rt, ql.lt);
    }
}; 

void export_point(){
    using namespace maths::geometry;

    class_<point, shared_ptr<point> >("Point", 
            init<const double&, const double&>())
        .def_pickle(PointPickler())
    ;
}


void export_quadrilateral(){
    using namespace maths::geometry;

    class_<quadrilateral, bases<Printable>, shared_ptr<quadrilateral> >("Quadrilateral", 
            init<const point&, const point&, const point&, 
            const point&, optional<const string&> >())
        .def_pickle(QuadrilateralPickler())
    ;
}

}
}

