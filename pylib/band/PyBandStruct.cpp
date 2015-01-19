/* 
 * File:   BandStruct.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 6, 2013, 9:25 AM
 */

#include "PyBandStruct.h"

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{

PyBandStruct::PyBandStruct(const Workers &workers, uint nn,  bool orthoBasis,
        bool calcEigV, const string &prefix): BandStruct(workers, nn,  orthoBasis, calcEigV, prefix)
{    
}

//void PyBandStruct::H(bp::object H, int ineigh){
//    shared_ptr<cxmat> HH = npy2mat<dcmplx>(H);
//    if (HH == 0){
//        throw runtime_error("In BandStruct::H(): cannot convert numpy array to C++ matrix");
//    }
//    BandStruct::H(HH, ineigh);
//}

void PyBandStruct::H(const cxmat& H, int ineigh){
    auto HH = std::make_shared<cxmat>(H);
    BandStruct::H(HH, ineigh);
}

//void PyBandStruct::S(bp::object S, int ineigh){
//    shared_ptr<cxmat> SS = npy2mat<dcmplx>(S);
//    if (SS == 0){
//        throw runtime_error("In BandStruct::HS(): cannot convert numpy array to C++ matrix");
//    }
//
//    BandStruct::S(SS, ineigh);
//}

void PyBandStruct::S(const cxmat& S, int ineigh){
    shared_ptr<cxmat> SS = make_shared<cxmat>(S);
    BandStruct::S(SS, ineigh);
}


void (PyBandStruct::*PyBandStruct_nb_set)(uint) = &PyBandStruct::nb;
uint (PyBandStruct::*PyBandStruct_nb_get)() = &PyBandStruct::nb;
void (PyBandStruct::*PyBandStruct_ne_set)(uint) = &PyBandStruct::ne;
uint (PyBandStruct::*PyBandStruct_ne_get)() = &PyBandStruct::ne;
void (PyBandStruct::*PyBandStruct_lv_set)(const lvec&) = &PyBandStruct::lv;
lvec (PyBandStruct::*PyBandStruct_lv_get)() = &PyBandStruct::lv;
void (PyBandStruct::*PyBandStruct_lc_1)(const lcoord&, int) = &PyBandStruct::lc;
void (PyBandStruct::*PyBandStruct_k_1)(const mat&) = &PyBandStruct::k;
//void (PyBandStruct::*PyBandStruct_H_1)(bp::object, int) = &PyBandStruct::H;
//void (PyBandStruct::*PyBandStruct_S_1)(bp::object, int) = &PyBandStruct::S;
void (PyBandStruct::*PyBandStruct_H_1)(const cxmat&, int) = &PyBandStruct::H;
void (PyBandStruct::*PyBandStruct_S_1)(const cxmat&, int) = &PyBandStruct::S;
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(PyBandStruct_save, save, 1, 2)
void export_BandStruct(){
    // ~~~~~~~~ To avoid nasty numpy segfault ~~~~~~~
    import_array();
    
    class_<BandStruct, bases<Printable>, shared_ptr<BandStruct>, noncopyable >("_BandStruct", 
            no_init)
    ;
    
    class_<PyBandStruct, bases<BandStruct>, shared_ptr<PyBandStruct>, noncopyable>("BandStruct",
            init<const Workers&, uint, optional<bool, bool, const string&> >())
        .add_property("nb", PyBandStruct_nb_get, PyBandStruct_nb_set)
        .add_property("ne", PyBandStruct_ne_get, PyBandStruct_ne_set)
        .add_property("lv", PyBandStruct_lv_get, PyBandStruct_lv_set)
        .def("lc", PyBandStruct_lc_1)
        .def("k", PyBandStruct_k_1)
        .def("H", PyBandStruct_H_1)
        .def("S", PyBandStruct_S_1)    
        .def("run", &BandStruct::run)
        .def("save", &BandStruct::save, PyBandStruct_save())
        .def("enableEigVec", &BandStruct::enableEigVec)
    ;
}

}
}


