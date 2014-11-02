/**
 * @file npyarma/npyarma.cpp
 * @date 
 * @author K M Masum Habib <masum.habib@gmail.com>.  
 *
 * @brief Automatic converters to-from numpy::array for armadillo. 
 * Adopted from automatic converters to-from python for blitz::Array for bob
 * project: http://www.idiap.ch/software/bob.
 *
 * Copyright (C) 2014 K M masum Habib.
 */

#include "npyarma/npyarma.h"



namespace qmicad{ namespace python{

namespace bp = boost::python;
using namespace utils::stds;
namespace ar = arma;
using std::complex;

template <> int ctype_to_npytype<bool>() { return NPY_BOOL; }

template <> int ctype_to_npytype<int8_t>() { return NPY_INT8; }
template <> int ctype_to_npytype<uint8_t>() { return NPY_UINT8; }
template <> int ctype_to_npytype<int16_t>() { return NPY_INT16; }
template <> int ctype_to_npytype<uint16_t>() { return NPY_UINT16; }
template <> int ctype_to_npytype<int32_t>() { return NPY_INT32; }
template <> int ctype_to_npytype<uint32_t>() { return NPY_UINT32; }
template <> int ctype_to_npytype<int64_t>() { return NPY_INT64; }
template <> int ctype_to_npytype<uint64_t>() { return NPY_UINT64; }
template <> int ctype_to_npytype<float>() { return NPY_FLOAT32; }
template <> int ctype_to_npytype<double>() { return NPY_FLOAT64; }
#ifdef NPY_FLOAT128
template <> int ctype_to_npytype<long double>() { return NPY_FLOAT128; }
#endif
template <> int ctype_to_npytype<std::complex<float> >() { return NPY_COMPLEX64; }
template <> int ctype_to_npytype<std::complex<double> >() { return NPY_COMPLEX128; }
#ifdef NPY_COMPLEX256
template <> int ctype_to_npytype<std::complex<long double> >() { return NPY_COMPLEX256; }
#endif

string type_to_string(int type){
    string type_name;
    switch(type){
        case NPY_BOOL: type_name = "NPY_BOOL"; break;
        case NPY_INT8: type_name = "NPY_INT8"; break;
        case NPY_UINT8: type_name = "NPY_UINT8"; break;
        case NPY_INT16: type_name = "NPY_INT16"; break;
        case NPY_UINT16: type_name = "NPY_UINT16"; break;
        case NPY_INT32: type_name = "NPY_INT32"; break;
        case NPY_UINT32: type_name = "NPY_UINT32"; break;
        case NPY_INT64: type_name = "NPY_INT64"; break;
        case NPY_UINT64: type_name = "NPY_UINT64"; break;
        case NPY_FLOAT32: type_name = "NPY_FLOAT32"; break;
        case NPY_FLOAT64: type_name = "NPY_FLOAT64"; break;
        #ifdef NPY_FLOAT128
        case NPY_FLOAT128: type_name = "NPY_FLOAT128"; break;
        #endif
        case NPY_COMPLEX64: type_name = "NPY_COMPLEX64"; break;
        case NPY_COMPLEX128: type_name = "NPY_COMPLEX128"; break;
        #ifdef NPY_COMPLEX256
        case NPY_COMPLEX256: type_name = "NPY_COMPLEX256"; break;
        #endif
        
        default: type_name = "Type not supported in C++."; break;
    }
    
    return type_name;
}

//bp::object construct(bp::object obj){
//    using namespace ar;
//    
//    cout << "got it" << endl;
//    
//    PyObject *objp = obj.ptr();
//    PyArrayObject* a = (PyArrayObject*)objp; //extract<PyArrayObject*>(obj);
//    
//    if (a == 0){
//        cout << "NULL" << endl;
//    }else{        
//        cout << "Not NULL" << endl;
//        cout << "DIM: " << a->nd << endl;
//        cout << "FLAGS: " << a->flags << endl;
//        cout << "Size: " << a->dimensions[0] << endl;
//        
//    }
//    
//    vec v((double*)a->data,  a->dimensions[0], false, true);
//    cout << "vec: " << endl << v << endl; 
//    v(0) = 5.0;
//    cout << "vec: " << endl << v << endl; 
//    
////    vec v2 = zeros<vec>(5);
////    v2(0) = 10;
////    cout << "vec2: " << endl << v2 << endl; 
////    npy_intp size[1] = {v2.n_elem};
////    PyObject * pyObj = PyArray_SimpleNewFromData(1, size, NPY_DOUBLE, (void*)v2.memptr());
//    npy_intp size[1] = {v.n_elem};
//    PyObject * pyObj = PyArray_SimpleNewFromData(1, size, NPY_DOUBLE, (void*)v.memptr());
////    PyObject * pyObj = (PyObject*)PyArray_SimpleNew(1, size, NPY_DOUBLE);
//    boost::python::handle<> handle( pyObj );
//    return bp::object(handle);
//    
//}

template <typename T>
void test_npy_to_mat(const ar::Mat<T> &value){
    cout << "In test_npy_to_mat(), value: " << endl << value << endl;
}

template <typename T>
ar::Mat<T> test_mat_to_npy(){
    ar::Mat<T> m(2,3);
    m << 1 << 2 << 3 << ar::endr
      << 4 << 5 << 6 << ar::endr;
    
    cout << "In test_mat_to_npy(), value: " << endl << m << endl;
    
    return m;
}

template <typename T>
void test_npy2mat(bp::object obj){
    shared_ptr<ar::Mat<T> > pmat = npy2mat<T>(obj);
    if (pmat !=0){
        cout << "In test_npy_to_mat(), value: " << endl << *pmat << endl;
    }
}

void export_npyarma(){
    import_array();

  /**
   * The following struct constructors will make sure we can input
   * arma::Col<T> in our bound C++ routines w/o needing to specify
   * special converters each time. The rvalue converters allow bp to
   * automatically map the following inputs:
   *
   * a) const arma::Col<T>& (pass by const reference)
   * b) arma::Col<T><T> (pass by value -- DO NEVER DO THIS!!!)
   *
   * Please note that the last case:
   * 
   * c) arma::Col<T>& (pass by non-const reference)
   *
   * is NOT covered by these converters. The reason being that because the
   * object may be changed, there is no way for bp to update the
   * original python object, in a sensible manner, at the return of the method.
   *
   * Avoid passing by non-const reference in your methods.
   */
    
   //col_from_npy<bool>();
   //col_from_npy<int8_t>();
   col_from_npy<int16_t>();
   col_from_npy<int32_t>();
   col_from_npy<int64_t>();
   col_from_npy<uint8_t>();
   col_from_npy<uint16_t>();
   col_from_npy<uint32_t>();
   col_from_npy<uint64_t>();
   col_from_npy<float>();
   col_from_npy<double>();
   //col_from_npy<long double>();
   col_from_npy<std::complex<float> >();
   col_from_npy<std::complex<double> >();
   //col_from_npy<std::complex<long double> >();

  
  /**
   * The following struct constructors will make C++ return values of type
   * arma::Col<T> to show up in the python side as numpy arrays.
   */
   //register_col_to_npy<bool>();
   //register_col_to_npy<int8_t>();
   register_col_to_npy<int16_t>();
   register_col_to_npy<int32_t>();
   register_col_to_npy<int64_t>();
   register_col_to_npy<uint8_t>();
   register_col_to_npy<uint16_t>();
   register_col_to_npy<uint32_t>();
   register_col_to_npy<uint64_t>();
   register_col_to_npy<float>();
   register_col_to_npy<double>();
   //register_col_to_npy<long double>();
   register_col_to_npy<std::complex<float> >();
   register_col_to_npy<std::complex<double> >();
   //register_col_to_npy<std::complex<long double> >();

   //row_from_npy<bool>();
   //row_from_npy<int8_t>();
   row_from_npy<int16_t>();
   row_from_npy<int32_t>();
   row_from_npy<int64_t>();
   row_from_npy<uint8_t>();
   row_from_npy<uint16_t>();
   row_from_npy<uint32_t>();
   row_from_npy<uint64_t>();
   row_from_npy<float>();
   row_from_npy<double>();
   //row_from_npy<long double>();
   row_from_npy<std::complex<float> >();
   row_from_npy<std::complex<double> >();
   //row_from_npy<std::complex<long double> >();

  
  /**
   * The following struct constructors will make C++ return values of type
   * arma::Col<T> to show up in the python side as numpy arrays.
   */
   //register_row_to_npy<bool>();
   //register_row_to_npy<int8_t>();
   register_row_to_npy<int16_t>();
   register_row_to_npy<int32_t>();
   register_row_to_npy<int64_t>();
   register_row_to_npy<uint8_t>();
   register_row_to_npy<uint16_t>();
   register_row_to_npy<uint32_t>();
   register_row_to_npy<uint64_t>();
   register_row_to_npy<float>();
   register_row_to_npy<double>();
   //register_row_to_npy<long double>();
   register_row_to_npy<std::complex<float> >();
   register_row_to_npy<std::complex<double> >();
   //register_row_to_npy<std::complex<long double> >();

   
   
   //mat_from_npy<bool>();
   //mat_from_npy<int8_t>();
   mat_from_npy<int16_t>();
   mat_from_npy<int32_t>();
   mat_from_npy<int64_t>();
   mat_from_npy<uint8_t>();
   mat_from_npy<uint16_t>();
   mat_from_npy<uint32_t>();
   mat_from_npy<uint64_t>();
   mat_from_npy<float>();
   mat_from_npy<double>();
   //mat_from_npy<long double>();
   mat_from_npy<std::complex<float> >();
   mat_from_npy<std::complex<double> >();
   //mat_from_npy<std::complex<long double> >();

  
  /**
   * The following struct constructors will make C++ return values of type
   * arma::Mat<T> to show up in the python side as numpy arrays.
   */
   //register_mat_to_npy<bool>();
   //register_mat_to_npy<int8_t>();
   register_mat_to_npy<int16_t>();
   register_mat_to_npy<int32_t>();
   register_mat_to_npy<int64_t>();
   register_mat_to_npy<uint8_t>();
   register_mat_to_npy<uint16_t>();
   register_mat_to_npy<uint32_t>();
   register_mat_to_npy<uint64_t>();
   register_mat_to_npy<float>();
   register_mat_to_npy<double>();
   //register_mat_to_npy<long double>();
   register_mat_to_npy<std::complex<float> >();
   register_mat_to_npy<std::complex<double> >();
   //register_mat_to_npy<std::complex<long double> >();
   
//   def("test_mat_to_npy", test_mat_to_npy<double>, " Tests arma::mat to numpy.array conversion.");
   def("test_mat_to_npy", test_mat_to_npy<complex<double> >, " Tests arma::mat to numpy.array conversion.");
//   def("test_npy_to_mat", test_npy_to_mat<double>, " Tests numpy.array to arma::mat conversion.");
   def("test_npy_to_mat", test_npy_to_mat<complex<double> >, " Tests numpy.array to arma::mat conversion.");
//   def("test_npy2mat", test_npy2mat<double>, " Tests numpy.array to arma::mat conversion.");
   def("test_npy2mat", test_npy2mat<complex<double> >, " Tests numpy.array to arma::mat conversion.");

   //   def("construct", construct, " Converts numpy.array to mat");


}


}}

