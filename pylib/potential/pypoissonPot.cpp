
/*
 * File:   pypoissonPot.cpp
 * Copyright (C) 2015  Mirza Elahi <mirza.monzur@gmail.com>
 *
 * Created on April 12, 2015, 2:56 AM
 */

#include "potential/poissonPot.h"
#include "boostpython.hpp"

/**
 * Python exporters.
 */
namespace qmicad {
namespace python {
using namespace potential;

/**
 * poisson potential
 */
//BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(poissonPot_VLR, VLR, 3, 5)

void export_poissonPot() {
    class_<poissonPot, bases<Potential>, shared_ptr < poissonPot >> ("poissonPot", init< vec, vec >())
            .def(init<double, double, double, double>())
            .def("setMaterialEps", &poissonPot::setMaterialEps)
            .def("setMaterialni", &poissonPot::setMaterialni)
            .def("setDoping", &poissonPot::setDoping)
            .def("setPotentialDirichlet", &poissonPot::setPotentialDirichlet)
            .def("setPhiByFroce", &poissonPot::setPhiByFroce)
            .def("generateGrad2LambdaMatrix", &poissonPot::generateGrad2LambdaMatrix)
            .def("setInitialGuess", &poissonPot::setInitialGuess)
            .def("calculateLambdaSingleIteration", &poissonPot::calculateLambdaSingleIteration)
            //.def("getRho", &poissonPot::getRho)
            .def("getPotentialSliceAlongZ", &poissonPot::getPotentialSliceAlongZ)
            .def("getPotentialSliceAlongX", &poissonPot::getPotentialSliceAlongX)
            .def("helperDoping", &poissonPot::helperDoping)
            
            
            //.enable_pickling()
            //        .def("addLinearRegion", &LinearPot::addLinearRegion)
            //        .add_property("NLR", &LinearPot::NLR)
            //        .def("VLR", &LinearPot::VLR, LinearPot_VLR()) 
            ;
}

}
}



