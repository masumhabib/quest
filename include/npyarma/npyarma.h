/**
 * @file python/npyarma.h
 * @date 
 * @author 
 *
 * @brief A
 *
 * Copyright (C) 
 */

#ifndef NPYARMA_H
#define NPYARMA_H

#include "python/boostpython.hpp"
#include "utils/std.hpp"
#include <numpy/arrayobject.h>
#include <armadillo>

namespace qmicad{ namespace python{

namespace bp = boost::python;
using namespace utils::stds;
namespace ar = arma;

template <typename T> int ctype_to_npytype() {
    throw runtime_error("Type not supported in C++.");
}


/**
 * Converter for boost::python::object & to Col<T>.
 */
template<typename T>
shared_ptr<ar::Col<T> > npy2col(bp::object obj){
    
    PyObject *objp = obj.ptr();
    PyArrayObject* a = (PyArrayObject*)objp;
    if(PyArray_Check(objp)){ // if this is a PyArray
        PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(objp);            
        if(ctype_to_npytype<T>() == PyArray_DESCR(arr)->type_num){ // if array has the same type
            if(arr->nd == 1){ // if array is 1 dimensional
                if(arr->flags &  NPY_F_CONTIGUOUS || arr->flags &  NPY_C_CONTIGUOUS){
                    
                    shared_ptr<ar::Col<T> > out(new ar::Col<T>((T*)a->data,  a->dimensions[0], false, true));
                    return out;
                }
            }
        }
    }    
}

/**
 * Converter for boost::python::object & to Row<T>.
 */
template<typename T>
shared_ptr<ar::Row<T> > npy2row(bp::object obj){
    
    PyObject *objp = obj.ptr();
    PyArrayObject* a = (PyArrayObject*)objp;
    if(PyArray_Check(objp)){ // if this is a PyArray
        PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(objp);            
        if(ctype_to_npytype<T>() == PyArray_DESCR(arr)->type_num){ // if array has the same type
            if(arr->nd == 1){ // if array is 1 dimensional
                if(arr->flags &  NPY_F_CONTIGUOUS || arr->flags &  NPY_C_CONTIGUOUS){
                    
                    shared_ptr<ar::Row<T> > out(new ar::Row<T>((T*)a->data,  a->dimensions[0], false, true));
                    return out;
                }
            }
        }
    }    
}


/**
 * Converter for boost::python::object & to Mat<T>.
 */
template<typename T>
shared_ptr<ar::Mat<T> > npy2mat(bp::object obj){

                    
    PyObject *objp = obj.ptr();
    PyArrayObject* a = (PyArrayObject*)objp;
    
    if(PyArray_Check(objp)){ // if this is a PyArray
        
        PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(objp);            
        
        if(ctype_to_npytype<T>() == PyArray_DESCR(arr)->type_num){ // if array has the same type
            if(arr->nd == 2){ 
                if(arr->flags &  NPY_F_CONTIGUOUS){
                    shared_ptr<ar::Mat<T> > out(new ar::Mat<T>((T*)a->data,  a->dimensions[0], a->dimensions[1], false, true));
                    return out;
                }else if (arr->flags &  NPY_C_CONTIGUOUS){
                    shared_ptr<ar::Mat<T> > out(new ar::Mat<T>(a->dimensions[0], a->dimensions[1]));
                    for (int i = 0; i < a->dimensions[0]; ++i){
                        for (int j = 0; j < a->dimensions[1]; ++j){
                            (*out)(i,j) = *static_cast<T*>(PyArray_GETPTR2(arr, i, j));
                        }
                    }        

                    return out;
                }
            }
        }
    }
    return  shared_ptr<ar::Mat<T> >(); // Failure, return null.
}



/**
 * Objects of this type create a binding between arma::Col<T> and
 * NumPy arrays. You can specify a NumPy array as a parameter to a
 * bound method that would normally receive a arma::Col<T> or a const
 * arma::Col<T>& and the conversion will just magically happen, as
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
     * Registers converter from numpy array into a arma::Col<T>
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

    if (obj_ptr == 0){
        throw runtime_error("Cannot convert numpy::array to arma::Col<T>. NULL object received.");
    }
      
    //black-magic required to setup the arma::Col<T> storage area
    void* storage = ((bp::converter::rvalue_from_python_storage<col_type>*)data)->storage.bytes;
    
    PyArrayObject *arr = reinterpret_cast<PyArrayObject*>(obj_ptr);
    
    //mounts the numpy memory at the "newly allocated" arma::Col<T>
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
  return (PyArrayObject*)PyArray_EMPTY(nd, dims, type, 1); //returns a fortran type ndarray.
}

/**
 * Objects of this type bind arma::Col<T> to numpy arrays. Your method
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

    //wrap new PyArray and then copy the data
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

/**
 * Objects of this type create a binding between arma::Row<T> and
 * NumPy arrays. You can specify a NumPy array as a parameter to a
 * bound method that would normally receive a arma::Row<T> or a const
 * arma::Row<T>& and the conversion will just magically happen, as
 * efficiently as possible.
 *
 * Please note that passing by value should be avoided as much as possible. In
 * this mode, the underlying method will still be able to alter the underlying
 * array storage area w/o being able to modify the array itself, causing a
 * gigantic mess. If you want to make something close to pass-by-value, just
 * pass by non-const reference instead.
 */
template <typename T> struct row_from_npy {
   
    typedef typename ar::Row<T> row_type;
    static const int N = 1;

    /**
     * Registers converter from numpy array into a arma::Row<T>
     */
    row_from_npy() {
      bp::converter::registry::push_back(&convertible, &construct, 
          bp::type_id<row_type>());
    }

    /**
     * This method will determine if the input python object is convertible into
     * a Row<T>
     */
    static void* convertible(PyObject* obj_ptr) {
        if(PyArray_Check(obj_ptr)){ // if this is a PyArray
            PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(obj_ptr);            
            if(ctype_to_npytype<T>() == PyArray_DESCR(arr)->type_num){ // if array has the same type
                if(arr->nd == N){ // if array is 1 dimensional
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

    if (obj_ptr == 0){
        throw runtime_error("Cannot convert numpy::array to arma::Row<T>. NULL object received.");
    }
      
    //black-magic required to setup the arma::Col<T> storage area
    void* storage = ((bp::converter::rvalue_from_python_storage<row_type>*)data)->storage.bytes;
    
    PyArrayObject *arr = reinterpret_cast<PyArrayObject*>(obj_ptr);
    
    //mounts the numpy memory at the "newly allocated" arma::Row<T>
    npy_intp shape[N];
    shape[0] = PyArray_DIMS(arr)[0];

    new (storage) row_type((T*)PyArray_DATA(arr), shape[0], false, true); //place operator
    data->convertible = storage;
}

};


/**
 * Objects of this type bind arma::Row<T> to numpy arrays. Your method
 * generates as output an object of this type and the object will be
 * automatically converted into a Numpy array.
 */
template <typename T> struct row_to_npy {
    
  static const int N = 1;
  typedef typename ar::Row<T> row_type;
  
  static PyObject* convert(const row_type& tv) {
    npy_intp dims[N];
    dims[0] = tv.n_elem;

    PyArrayObject* retval = make_pyarray(N, dims, ctype_to_npytype<T>());

    //wrap new PyArray and then copy the data
    row_type rowdest((T*)PyArray_DATA(retval), dims[0], false, true);
    rowdest = tv;

    return reinterpret_cast<PyObject*>(retval);
  }

  static const PyTypeObject* get_pytype() { return &PyArray_Type; }

};

template <typename T>
void register_row_to_npy() {
  bp::to_python_converter<typename ar::Row<T>, row_to_npy<T>
#if defined BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                          ,true
#endif
              >();
}

/**
 * Objects of this type create a binding between arma::Mat<T> and
 * NumPy arrays. You can specify a NumPy array as a parameter to a
 * bound method that would normally receive a arma::Mat<T> or a const
 * arma::Mat<T>& and the conversion will just magically happen, as
 * efficiently as possible.
 *
 * Please note that passing by value should be avoided as much as possible. In
 * this mode, the underlying method will still be able to alter the underlying
 * array storage area w/o being able to modify the array itself, causing a
 * gigantic mess. If you want to make something close to pass-by-value, just
 * pass by non-const reference instead.
 */
template <typename T> struct mat_from_npy {
   
    typedef typename ar::Mat<T> mat_type;
    static const int N = 2;

    /**
     * Registers converter from numpy array into a arma::Mat<T>
     */
    mat_from_npy() {
      bp::converter::registry::push_back(&convertible, &construct, 
          bp::type_id<mat_type>());
    }

    /**
     * This method will determine if the input python object is convertible into
     * a Mat<T>
     */
    static void* convertible(PyObject* obj_ptr) {
        if(PyArray_Check(obj_ptr)){ // if this is a PyArray
            PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(obj_ptr);            
            if(ctype_to_npytype<T>() == PyArray_DESCR(arr)->type_num){ // if array has the same type
                if(arr->nd == N){ // if array is 2 dimensional
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

    if (obj_ptr == 0){
        throw runtime_error("Cannot convert numpy::array to arma::Mat<T>. NULL object received.");
    }
      
    //black-magic required to setup the arma::Col<T> storage area
    void* storage = ((bp::converter::rvalue_from_python_storage<mat_type>*)data)->storage.bytes;
    
    PyArrayObject *arr = reinterpret_cast<PyArrayObject*>(obj_ptr);
    
    //mounts the numpy memory at the "newly allocated" arma::Mat<T>
    npy_intp shape[N], stride[N];
    shape[0] = PyArray_DIMS(arr)[0];
    shape[1] = PyArray_DIMS(arr)[1];
    
    stride[0] = (PyArray_STRIDES(arr)[0]/sizeof(T));
    stride[1] = (PyArray_STRIDES(arr)[1]/sizeof(T));

    if (arr->flags &  NPY_F_CONTIGUOUS){
        // We have the right array type. So just create a mat_type
        new (storage) mat_type((T*)PyArray_DATA(arr), shape[0], shape[1], false, true); //place operator
    }else {
        // Array is C contiguous. We need to copy.
        new (storage) mat_type(shape[0], shape[1]); //place operator
        mat_type *mp = (mat_type*)storage;
        for (int i = 0; i < shape[0]; ++i){
            for (int j = 0; j < shape[1]; ++j){
                (*mp)(i,j) = *static_cast<T*>(PyArray_GETPTR2(arr, i, j));
            }
        }        
    }
    
    data->convertible = storage;
}

};


/**
 * Objects of this type bind arma::Mat<T> to numpy arrays. Your method
 * generates as output an object of this type and the object will be
 * automatically converted into a Numpy array.
 */
template <typename T> struct mat_to_npy {
    
  static const int N = 2;
  typedef typename ar::Mat<T> mat_type;
  
  static PyObject* convert(const mat_type& tv) {
    npy_intp dims[N];
    dims[0] = tv.n_rows;
    dims[1] = tv.n_cols;

    PyArrayObject* retval = make_pyarray(N, dims, ctype_to_npytype<T>());

    //wrap new PyArray and then copy the data
    mat_type matdest((T*)PyArray_DATA(retval), dims[0], dims[1], false, true);
    matdest = tv;

    return reinterpret_cast<PyObject*>(retval);
  }

  static const PyTypeObject* get_pytype() { return &PyArray_Type; }

};

template <typename T>
void register_mat_to_npy() {
  bp::to_python_converter<typename ar::Mat<T>, mat_to_npy<T>
#if defined BOOST_PYTHON_SUPPORTS_PY_SIGNATURES
                          ,true
#endif
              >();
}


}}

#endif /* NPYARMA_H */
