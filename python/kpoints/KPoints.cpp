/* 
 * File:   KPoints.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 17, 2014, 12:21 AM
 */

#include "kpoints/KPoints.h"
#include "python/boostpython.hpp"


/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace kpoints;

void export_KPoints(){
    void (KPoints::*KPoints_addKLine1)(const point&, const point&, double) = &KPoints::addKLine;
    void (KPoints::*KPoints_addKLine2)(const point&, const point&, uint) = &KPoints::addKLine;
    void (KPoints::*KPoints_addKRect1)(const point&, const point&, double, double) = &KPoints::addKRect;
    void (KPoints::*KPoints_addKRect2)(const point&, const point&, uint, uint) = &KPoints::addKRect;
    
    class_<KPoints, bases<Printable>, shared_ptr<KPoints>, noncopyable>("KPoints",
            init<optional<const string&> >()) 
        .def("addKPoint", &KPoints::addKPoint)
        .def("addKLine", KPoints_addKLine1, "Creates k-grid along a line for a given interval.")
        .def("addKLine", KPoints_addKLine2, "Creates k-grid along a line for a given number of k-points.")
        .def("addKRect", KPoints_addKRect1, "Creates a rectangular k-grid for a given interval.")
        .def("addKRect", KPoints_addKRect2, "Creates a rectangular k-grid for a given number of k-points.")
        .add_property("kp", &KPoints::kp, "k-points property, readonly.")
        .add_property("N", &KPoints::N, "Total number of k-points, readonly.")
    ;
}
}
}



