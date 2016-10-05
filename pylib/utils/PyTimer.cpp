/* 
 * File:   Timer.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 11, 2014, 12:43 AM
 */

#include "boostpython.hpp"
#include "utils/Timer.h"
#include "maths/geometry.hpp"


namespace quest{
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
