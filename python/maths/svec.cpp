/* 
 * File:   svec.cpp
 * Copyright (C) 2013  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 7, 2013, 9:51 AM
 * 
 * Description: Spetial vector type definition.
 * 
 */

#include "maths/svec.h"
#include "python/boostpython.hpp"


namespace qmicad{
namespace python{


void export_svec(){
    using namespace maths::spvec;
    
    /**
     * Spatial vector/position vector.
     */
    //class_<svec, shared_ptr<svec> >("svec")
    //;
}

void export_pvec(){
    using namespace maths::spvec;
    
    /**
     * Position vector. Just a wrapper of svec.
     */
    class_<SVec, /*bases<svec>,*/ shared_ptr<SVec> >("SVec", 
            init<optional<double, double, double> >())
        .def(init<const svec&>())    
        .add_property("X", &SVec::getx, &SVec::setx)
        .add_property("Y", &SVec::gety, &SVec::sety)
        .add_property("Z", &SVec::getz, &SVec::setz)
    ;    
}

}
}



