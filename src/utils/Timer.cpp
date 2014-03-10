/* 
 * File:   Timer.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 * 
 * Created on February 11, 2014, 12:43 AM
 */

#include "Timer.h"
#include "../maths/geometry.hpp"
#include "../python/boostpython.hpp"


namespace qmicad{
namespace python{

/**
 * Wall clock python exporter.
 */
void export_Timer(){
    using namespace maths::geometry;
    
    class_<Timer, bases<Printable>, shared_ptr<Timer> >("Timer", 
            init<optional<const string&> >())
        .def("tic", &Timer::tic)
        .def("toc", &Timer::toc)
    ;

}

}
}