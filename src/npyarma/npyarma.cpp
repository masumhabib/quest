/**
 * @file python/blitz_numpy.cc
 * @date Mon Sep 26 11:47:30 2011 +0200
 * @author Andre Anjos <andre.anjos@idiap.ch>
 *
 * @brief Automatic converters to-from python for blitz::Array's
 *
 * Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland
 */

#include "npyarma/npyarma.h"
#include "utils/std.hpp"



namespace qmicad{ namespace python{

namespace bp = boost::python;
using namespace utils::stds;
namespace ar = arma;

template <typename T> int ctype_to_npytype() {
    throw runtime_error("Type not supported in C++");
    // throw error
    //PYTHON_ERROR(TypeError, "unsupported C/C++ type (%s)", stringize<T>());
}

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

//template<typename T, int N>
//void npy_copy_cast(blitz::Array<T,N>& bz, PyArrayObject* arrobj) {
//  PYTHON_ERROR(TypeError, "unsupported number of dimensions: %d", N);
//}

//template<typename T>
//static void npy_copy_cast(blitz::Array<T,1>& bz, PyArrayObject* arrobj) {
//  for (int i=0; i<PyArray_DIM(arrobj,0); ++i)
//    bz(i) = *static_cast<T*>(PyArray_GETPTR1(arrobj, i));
//}

//template<typename T>
//static void npy_copy_cast(blitz::Array<T,2>& bz, PyArrayObject* arrobj) {
//  for (int i=0; i<PyArray_DIM(arrobj,0); ++i)
//    for (int j=0; j<PyArray_DIM(arrobj,1); ++j)
//      bz(i,j) = *static_cast<T*>(PyArray_GETPTR2(arrobj, i, j));
//}

//template<typename T>
//static void npy_copy_cast(blitz::Array<T,3>& bz, PyArrayObject* arrobj) {
//  for (int i=0; i<PyArray_DIM(arrobj,0); ++i)
//    for (int j=0; j<PyArray_DIM(arrobj,1); ++j)
//      for (int k=0; k<PyArray_DIM(arrobj,2); ++k)
//        bz(i,j,k) = *static_cast<T*>(PyArray_GETPTR3(arrobj, i, j, k));
//}

//template<typename T>
//static void npy_copy_cast(blitz::Array<T,4>& bz, PyArrayObject* arrobj) {
//  for (int i=0; i<PyArray_DIM(arrobj,0); ++i)
//    for (int j=0; j<PyArray_DIM(arrobj,1); ++j)
//      for (int k=0; k<PyArray_DIM(arrobj,2); ++k)
//        for (int l=0; l<PyArray_DIM(arrobj,3); ++l)
//          bz(i,j,k,l) = *static_cast<T*>(PyArray_GETPTR4(arrobj, i, j, k, l));
//}

/**
 * Objects of this type create a binding between blitz::Array<T,N> and
 * NumPy arrays. You can specify a NumPy array as a parameter to a
 * bound method that would normally receive a blitz::Array<T,N> or a const
 * blitz::Array<T,N>& and the conversion will just magically happen, as
 * efficiently as possible.
 *
 * Please note that passing by value should be avoided as much as possible. In
 * this mode, the underlying method will still be able to alter the underlying
 * array storage area w/o being able to modify the array itself, causing a
 * gigantic mess. If you want to make something close to pass-by-value, just
 * pass by non-const reference instead.
 */
template <typename T> struct col_from_npy {
   
    typedef typename ar::Col<T> col_type;
    static const int N = 1;

    /**
     * Registers converter from numpy array into a ar::Col<T>
     */
    col_from_npy() {
      bp::converter::registry::push_back(&convertible, &construct, 
          bp::type_id<col_type>());
    }

    /**
     * This method will determine if the input python object is convertible into
     * a Col<T>
     */
    static void* convertible(PyObject* obj_ptr) {
        PyObject *out = 0;
        if(PyArray_Check(obj_ptr)){ // if this is a PyArray
            PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(obj_ptr);            
            if(ctype_to_npytype<T>() == PyArray_DESCR(arr)->type_num){ // if array has the same type
                if(arr->nd == 1){ // if array is 1 dimensional
                    if(arr->flags &  NPY_F_CONTIGUOUS || arr->flags &  NPY_C_CONTIGUOUS){
                        return obj_ptr;
                    }
                }
            }
        }

        return 0;
    }

  /**
   * This method will finally construct the C++ element out of the python
   * object that was input. Please note that when bp reaches this
   * method, the object has already been checked for convertibility.
   */
static void construct(PyObject* obj_ptr,
        bp::converter::rvalue_from_python_stage1_data* data) {
    assert(obj_ptr);
    
    if (obj_ptr == 0){
        cout << "ERROR" << endl;
    }
      
    //black-magic required to setup the blitz::Array<> storage area
    void* storage = ((bp::converter::rvalue_from_python_storage<col_type>*)data)->storage.bytes;
    
    PyArrayObject *arr = reinterpret_cast<PyArrayObject*>(obj_ptr);
    
    //mounts the numpy memory at the "newly allocated" blitz::Array
    npy_intp shape[N];
    shape[0] = PyArray_DIMS(arr)[0];

    new (storage) col_type((T*)PyArray_DATA(arr), shape[0], false, true); //place operator
    data->convertible = storage;
}

};

/**
 * Avoids the big number of warnings...
 */
static PyArrayObject* make_pyarray(int nd, npy_intp* dims, int type) {
  return (PyArrayObject*)PyArray_SimpleNew(nd, dims, type);
}

/**
 * Objects of this type bind blitz::Array<T,N> to numpy arrays. Your method
 * generates as output an object of this type and the object will be
 * automatically converted into a Numpy array.
 */
template <typename T> struct col_to_npy {
    
  static const int N = 1;
  typedef typename ar::Col<T> col_type;
  
  static PyObject* convert(const col_type& tv) {
    npy_intp dims[N];
    dims[0] = tv.n_elem;

    PyArrayObject* retval = make_pyarray(N, dims, ctype_to_npytype<T>());

    //wrap new PyArray in a blitz layer and then copy the data
    col_type coldest((T*)PyArray_DATA(retval), dims[0], false, true);
    coldest = tv;

    return reinterpret_cast<PyObject*>(retval);
  }

  static const PyTypeObject* get_pytype() { return &PyArray_Type; }

};

template <typename T>
void register_col_to_npy() {
  bp::to_python_converter<typename ar::Col<T>, col_to_npy<T>
#if defined BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                          ,true
#endif
              >();
}


bp::object construct(bp::object obj){
    using namespace ar;
    
    cout << "got it" << endl;
    
    PyObject *objp = obj.ptr();
    PyArrayObject* a = (PyArrayObject*)objp; //extract<PyArrayObject*>(obj);
    
    if (a == 0){
        cout << "NULL" << endl;
    }else{        
        cout << "Not NULL" << endl;
        cout << "DIM: " << a->nd << endl;
        cout << "FLAGS: " << a->flags << endl;
        cout << "Size: " << a->dimensions[0] << endl;
        
    }
    
    vec v((double*)a->data,  a->dimensions[0], false, true);
    cout << "vec: " << endl << v << endl; 
    v(0) = 5.0;
    cout << "vec: " << endl << v << endl; 
    
//    vec v2 = zeros<vec>(5);
//    v2(0) = 10;
//    cout << "vec2: " << endl << v2 << endl; 
//    npy_intp size[1] = {v2.n_elem};
//    PyObject * pyObj = PyArray_SimpleNewFromData(1, size, NPY_DOUBLE, (void*)v2.memptr());
    npy_intp size[1] = {v.n_elem};
    PyObject * pyObj = PyArray_SimpleNewFromData(1, size, NPY_DOUBLE, (void*)v.memptr());
//    PyObject * pyObj = (PyObject*)PyArray_SimpleNew(1, size, NPY_DOUBLE);
    boost::python::handle<> handle( pyObj );
    return bp::object(handle);
    
}

void test(/*const*/ ar::vec &v){
    cout << "vec: " << endl << v;
    
    v(0) = 10;
    
    cout << "vec: " << endl << v;
}


void export_npyarma(){
    import_array();
    def("construct", construct, " Converts numpy.array to mat");

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
   
   def("test", test, " Converts numpy.array to mat");

}


//void bind_core_bz_numpy () {
//  /**
//   * The following struct constructors will make sure we can input
//   * blitz::Array<T,N> in our bound C++ routines w/o needing to specify
//   * special converters each time. The rvalue converters allow bp to
//   * automatically map the following inputs:
//   *
//   * a) const blitz::Array<T,N>& (pass by const reference)
//   * b) blitz::Array<T,N> (pass by value -- DO NEVER DO THIS!!!)
//   *
//   * Please note that the last case:
//   * 
//   * c) blitz::Array<T,N>& (pass by non-const reference)
//   *
//   * is NOT covered by these converters. The reason being that because the
//   * object may be changed, there is no way for bp to update the
//   * original python object, in a sensible manner, at the return of the method.
//   *
//   * Avoid passing by non-const reference in your methods.
//   */
//#  define BOOST_PP_LOCAL_LIMITS (1, BOB_MAX_DIM)
//#  define BOOST_PP_LOCAL_MACRO(D) \
//   bz_from_npy<bool,D>();\
//   bz_from_npy<int8_t,D>();\
//   bz_from_npy<int16_t,D>();\
//   bz_from_npy<int32_t,D>();\
//   bz_from_npy<int64_t,D>();\
//   bz_from_npy<uint8_t,D>();\
//   bz_from_npy<uint16_t,D>();\
//   bz_from_npy<uint32_t,D>();\
//   bz_from_npy<uint64_t,D>();\
//   bz_from_npy<float,D>();\
//   bz_from_npy<double,D>();\
//   bz_from_npy<long double,D>();\
//   bz_from_npy<std::complex<float>,D>();\
//   bz_from_npy<std::complex<double>,D>();\
//   bz_from_npy<std::complex<long double>,D>();
//#  include BOOST_PP_LOCAL_ITERATE()
//  
//  /**
//   * The following struct constructors will make C++ return values of type
//   * blitz::Array<T,N> to show up in the python side as numpy arrays.
//   */
//#  define BOOST_PP_LOCAL_LIMITS (1, BOB_MAX_DIM)
//#  define BOOST_PP_LOCAL_MACRO(D) \
//   register_bz_to_npy<bool,D>();\
//   register_bz_to_npy<int8_t,D>();\
//   register_bz_to_npy<int16_t,D>();\
//   register_bz_to_npy<int32_t,D>();\
//   register_bz_to_npy<int64_t,D>();\
//   register_bz_to_npy<uint8_t,D>();\
//   register_bz_to_npy<uint16_t,D>();\
//   register_bz_to_npy<uint32_t,D>();\
//   register_bz_to_npy<uint64_t,D>();\
//   register_bz_to_npy<float,D>();\
//   register_bz_to_npy<double,D>();\
//   register_bz_to_npy<long double,D>();\
//   register_bz_to_npy<std::complex<float>,D>();\
//   register_bz_to_npy<std::complex<double>,D>();\
//   register_bz_to_npy<std::complex<long double>,D>();
//#  include BOOST_PP_LOCAL_ITERATE()
//}



}}

