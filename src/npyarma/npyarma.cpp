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
//#include "include/ndarray.h"


namespace qmicad{ namespace python{

namespace bp = boost::python;
using namespace utils::stds;
using namespace arma;

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
//template <typename T> struct mat_from_npy {
//   
//    typedef typename ar::Mat<T> mat_type;
//    static const N = 2;
//
//  /**
//   * Registers converter from numpy array into a blitz::Array<T,N>
//   */
//  mat_from_npy() {
//    bp::converter::registry::push_back(&convertible, &construct, 
//        bp::type_id<mat_type>());
//  }
//
//  /**
//   * This method will determine if the input python object is convertible into
//   * a Mat<T>
//   */
//  static void* convertible(PyObject* obj_ptr) {
//    bp::handle<> hdl(bp::borrowed(bp::allow_null(obj_ptr)));
//    bp::object obj(hdl);
//
//    typeinfo tinfo(getElementType<T>(), N);
//
//    convert_t result = convertible_to(obj, tinfo, false, true);
//
//    // we cannot afford copying, only referencing.
//    if (result == BYREFERENCE) return obj_ptr;
//
//    // but, if the user passed an array of the right type, but we still need to
//    // copy, warn the user as this is a tricky case to debug.
//    PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(obj_ptr);
//    if (result == WITHARRAYCOPY && 
//        ctype_to_num<T>() == PyArray_DESCR(arr)->type_num) {
//      PYTHON_ERROR(RuntimeError, "The bindings you are trying to use to this C++ method require a numpy.ndarray -> blitz::Array<%s,%d> conversion, but the array you passed, despite the correct type, is not C-style contiguous and/or properly aligned, so I cannot automatically wrap it. You can check this by yourself by printing the flags on such a variable with the command 'print(<varname>.flags)'. The only way to circumvent this problem, from python, is to create a copy the variable by issuing '<varname>.copy()' before calling the bound method. Otherwise, if you wish the copy to be executed automatically, you have to re-bind the method to use our custom 'const_ndarray' type.", stringize<T>(), N);
//    }
//
//    return 0;
//  }
//
//  /**
//   * This method will finally construct the C++ element out of the python
//   * object that was input. Please note that when bp reaches this
//   * method, the object has already been checked for convertibility.
//   */
//  static void construct(PyObject* obj_ptr,
//      bp::converter::rvalue_from_python_stage1_data* data) {
//
//    //black-magic required to setup the blitz::Array<> storage area
//    void* storage = ((bp::converter::rvalue_from_python_storage<mat_type>*)data)->storage.bytes;
//
//    PyArrayObject *arr = reinterpret_cast<PyArrayObject*>(obj_ptr);
//    
//    //mounts the numpy memory at the "newly allocated" blitz::Array
//    shape_type shape;
//    shape_type stride;
//    
//    for (int k=0; k<N; ++k) {
//      shape[k] = PyArray_DIMS(arr)[k];
//      stride[k] = (PyArray_STRIDES(arr)[k]/sizeof(T));
//    }
//    new (storage) mat_type((T*)PyArray_DATA(arr), shape[0], shape[1], false, true); //place operator
//    data->convertible = storage;
//
//  }
//
//};
//
/////**
//// * Avoids the big number of warnings...
//// */
////static PyArrayObject* make_pyarray(int nd, npy_intp* dims, int type) {
////  return (PyArrayObject*)PyArray_SimpleNew(nd, dims, type);
////}
////
///**
// * Objects of this type bind blitz::Array<T,N> to numpy arrays. Your method
// * generates as output an object of this type and the object will be
// * automatically converted into a Numpy array.
// */
//template <typename T> struct mat_to_npy {
//
//    
//  typedef typename ar::Mat<T> mat_type;
//  typedef typename npy_intp shape_type[N];
//
//  static const int N = 2;
//  
//  static PyObject* convert(const mat_type& tv) {
//    npy_intp dims[N];
//    for (int i=0; i<N; ++i) dims[i] = tv.extent(i);
//
//    PyArrayObject* retval = make_pyarray(N, dims, ctype_to_num<T>());
//
//    //wrap new PyArray in a blitz layer and then copy the data
//    shape_type shape=0;
//    for (int k=0; k<PyArray_NDIM(retval); ++k) shape[k] = PyArray_DIMS(retval)[k];
//    shape_type stride=0;
//    for (int k=0; k<PyArray_NDIM(retval); ++k) stride[k] = (PyArray_STRIDES(retval)[k]/sizeof(T));
//    mat_type bzdest((T*)PyArray_DATA(retval), shape, stride, blitz::neverDeleteData);
//    bzdest = tv;
//
//    return reinterpret_cast<PyObject*>(retval);
//  }
//
//  static const PyTypeObject* get_pytype() { return &PyArray_Type; }
//
//};

//template <typename T, int N>
//void register_bz_to_npy() {
//  bp::to_python_converter<typename blitz::Array<T,N>, bz_to_npy<T,N>
//#if defined BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
//                          ,true
//#endif
//              >();
//}
//
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


bp::object construct(bp::object obj){

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
    
    vec v2 = zeros<vec>(5);
    v2(0) = 10;
    cout << "vec2: " << endl << v2 << endl; 
//     Py_Initialize();
    npy_intp size[1] = {v2.n_elem};
//    int *data = new int[3];
    PyObject * pyObj = PyArray_SimpleNewFromData(1, size, NPY_DOUBLE, (void*)v2.memptr());
//    PyObject * pyObj = (PyObject*)PyArray_SimpleNew(1, size, NPY_DOUBLE);
    boost::python::handle<> handle( pyObj );
    //int n;
    //cout << "Hey";
    //cin >> n; 
    return bp::object(handle);
    //return bp::object();
//    cout  << obj.attr("dtype") << endl;
//  //black-magic required to setup the blitz::Array<> storage area
//  void* storage = ((bp::converter::rvalue_from_python_storage<mat_type>*)data)->storage.bytes;
//
//  PyArrayObject *arr = reinterpret_cast<PyArrayObject*>(obj_ptr);
//
//  //mounts the numpy memory at the "newly allocated" blitz::Array
//  shape_type shape;
//  shape_type stride;
//
//  for (int k=0; k<N; ++k) {
//    shape[k] = PyArray_DIMS(arr)[k];
//    stride[k] = (PyArray_STRIDES(arr)[k]/sizeof(T));
//  }
//  new (storage) mat_type((T*)PyArray_DATA(arr), shape[0], shape[1], false, true); //place operator
//  data->convertible = storage;

}


void export_npyarma(){
    import_array();
    def("construct", construct, " Converts numpy.array to mat");
}

}}

