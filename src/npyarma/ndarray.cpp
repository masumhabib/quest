/**
 * @file python/ndarray.cc
 * @date Tue Jan 18 17:07:26 2011 +0100
 * @author André Anjos <andre.anjos@idiap.ch>
 *
 * @brief Implementation of the ndarray class
 *
 * Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland
 */

//#include <boost/python/numeric.hpp>
//#include <boost/format.hpp>
//#include <stdexcept>
//#include <dlfcn.h>
//
//#include "ndarray.h"
//
//
//#define TP_ARRAY(x) ((PyArrayObject*)x.ptr())
//#define TP_OBJECT(x) (x.ptr())
//
//#define NUMPY17_API 0x00000007
//#define NUMPY16_API 0x00000006
//#define NUMPY14_API 0x00000004
//
//#if PY_VERSION_HEX >= 0x03000000
//static void* wrap_import_array() {
//  import_array();
//  return 0;
//}
//#else
//static void wrap_import_array() {
//  import_array();
//  return;
//}
//#endif
//
//namespace qmicad{ namespace python{
//
//void setup_python(const char* module_docstring) {
//
//  // Documentation options
//  if (module_docstring) bp::scope().attr("__doc__") = module_docstring;
//
//  // Required for logging C++ <-> Python interaction
//  if (!PyEval_ThreadsInitialized()) PyEval_InitThreads();
//
//  // Gets the current dlopenflags and save it
//  PyThreadState* tstate = PyThreadState_Get();
//  int old_value = tstate->interp->dlopenflags;
//
//  // Unsets the RTLD_GLOBAL flag
//  tstate->interp->dlopenflags = old_value & (~RTLD_GLOBAL);
//
//  // Loads numpy with the RTLD_GLOBAL flag unset
//  wrap_import_array();
//
//  // Resets the RTLD_GLOBAL flag
//  tstate->interp->dlopenflags = old_value;
//
//  //Sets the bp::numeric::array interface to use numpy.ndarray
//  //as basis. This is not strictly required, but good to set as a baseline.
//  bp::numeric::array::set_module_and_type("numpy", "ndarray");
//
//  // Make sure we are not running against the wrong version of NumPy
//  if (NPY_VERSION != PyArray_GetNDArrayCVersion()) {
//    PYTHON_ERROR(ImportError, "module compiled against ABI version 0x%08x but this version of numpy is 0x%08x - make sure you compile and execute against the same or compatible versions", (int) NPY_VERSION, (int) PyArray_GetNDArrayCVersion());
//  }
//
//#if NPY_FEATURE_VERSION >= NUMPY14_API /* NumPy C-API version >= 1.4 */
//  if (NPY_FEATURE_VERSION > PyArray_GetNDArrayCFeatureVersion()) {
//    PYTHON_ERROR(ImportError, "module compiled against API version 0x%08x but this version of numpy is 0x%08x - make sure you compile and execute against the same or compatible versions", (int) NPY_FEATURE_VERSION, (int) PyArray_GetNDArrayCFeatureVersion());
//  }
//#endif
//
//}
//
///***************************************************************************
// * Dtype (PyArray_Descr) manipulations                                     *
// ***************************************************************************/
//
//int type_to_num(ElementType type) {
//
//  switch(type) {
//    case t_bool:
//      return NPY_BOOL;
//    case t_int8:
//      return NPY_INT8;
//    case t_int16:
//      return NPY_INT16;
//    case t_int32:
//      return NPY_INT32;
//    case t_int64:
//      return NPY_INT64;
//    case t_uint8:
//      return NPY_UINT8;
//    case t_uint16:
//      return NPY_UINT16;
//    case t_uint32:
//      return NPY_UINT32;
//    case t_uint64:
//      return NPY_UINT64;
//    case t_float32:
//      return NPY_FLOAT32;
//    case t_float64:
//      return NPY_FLOAT64;
//#ifdef NPY_FLOAT128
//    case t_float128:
//      return NPY_FLOAT128;
//#endif
//    case t_complex64:
//      return NPY_COMPLEX64;
//    case t_complex128:
//      return NPY_COMPLEX128;
//#ifdef NPY_COMPLEX256
//    case t_complex256:
//      return NPY_COMPLEX256;
//#endif
//    default:
//      PYTHON_ERROR(TypeError, "unsupported C++ element type (%s)", stringize(type));
//  }
//
//}
//
//static ElementType signed_integer_type(int bits) {
//  switch(bits) {
//    case 8:
//      return t_int8;
//    case 16:
//      return t_int16;
//    case 32:
//      return t_int32;
//    case 64:
//      return t_int64;
//    default:
//      PYTHON_ERROR(TypeError, "unsupported signed integer element type with %d bits", bits);
//  }
//}
//
//static ElementType unsigned_integer_type(int bits) {
//  switch(bits) {
//    case 8:
//      return t_uint8;
//    case 16:
//      return t_uint16;
//    case 32:
//      return t_uint32;
//    case 64:
//      return t_uint64;
//    default:
//      PYTHON_ERROR(TypeError, "unsupported unsigned integer element type with %d bits", bits);
//  }
//}
//
//ElementType num_to_type(int num) {
//  switch(num) {
//    case NPY_BOOL:
//      return t_bool;
//
//    //signed integers
//    case NPY_BYTE:
//      return signed_integer_type(NPY_BITSOF_CHAR);
//    case NPY_SHORT:
//      return signed_integer_type(NPY_BITSOF_SHORT);
//    case NPY_INT:
//      return signed_integer_type(NPY_BITSOF_INT);
//    case NPY_LONG:
//      return signed_integer_type(NPY_BITSOF_LONG);
//    case NPY_LONGLONG:
//      return signed_integer_type(NPY_BITSOF_LONGLONG);
//
//    //unsigned integers
//    case NPY_UBYTE:
//      return unsigned_integer_type(NPY_BITSOF_CHAR);
//    case NPY_USHORT:
//      return unsigned_integer_type(NPY_BITSOF_SHORT);
//    case NPY_UINT:
//      return unsigned_integer_type(NPY_BITSOF_INT);
//    case NPY_ULONG:
//      return unsigned_integer_type(NPY_BITSOF_LONG);
//    case NPY_ULONGLONG:
//      return unsigned_integer_type(NPY_BITSOF_LONGLONG);
//
//    //floats
//    case NPY_FLOAT32:
//      return t_float32;
//    case NPY_FLOAT64:
//      return t_float64;
//#ifdef NPY_FLOAT128
//    case NPY_FLOAT128:
//      return t_float128;
//#endif
//
//    //complex
//    case NPY_COMPLEX64:
//      return t_complex64;
//    case NPY_COMPLEX128:
//      return t_complex128;
//#ifdef NPY_COMPLEX256
//    case NPY_COMPLEX256:
//      return t_complex256;
//#endif
//
//    default:
//      PYTHON_ERROR(TypeError, "unsupported NumPy element type (%d)", num);
//  }
//
//}
//
//template <> int ctype_to_num<bool>(void)
//{ return NPY_BOOL; }
//
////! @cond SKIP_DOXYGEN_WARNINGS
//template <> int ctype_to_num<int8_t>(void)
//{ return NPY_INT8; }
//template <> int ctype_to_num<uint8_t>(void)
//{ return NPY_UINT8; }
//template <> int ctype_to_num<int16_t>(void)
//{ return NPY_INT16; }
//template <> int ctype_to_num<uint16_t>(void)
//{ return NPY_UINT16; }
//template <> int ctype_to_num<int32_t>(void)
//{ return NPY_INT32; }
//template <> int ctype_to_num<uint32_t>(void)
//{ return NPY_UINT32; }
//template <> int ctype_to_num<int64_t>(void)
//{ return NPY_INT64; }
//template <> int ctype_to_num<uint64_t>(void)
//{ return NPY_UINT64; }
//template <> int ctype_to_num<float>(void)
//{ return NPY_FLOAT32; }
//template <> int ctype_to_num<double>(void)
//{ return NPY_FLOAT64; }
//#ifdef NPY_FLOAT128
//template <> int ctype_to_num<long double>(void)
//{ return NPY_FLOAT128; }
//#endif
//template <> int ctype_to_num<std::complex<float> >(void)
//{ return NPY_COMPLEX64; }
//template <> int ctype_to_num<std::complex<double> >(void)
//{ return NPY_COMPLEX128; }
//#ifdef NPY_COMPLEX256
//template <> int ctype_to_num<std::complex<long double> >(void)
//{ return NPY_COMPLEX256; }
//#endif
////! @endcond SKIP_DOXYGEN_WARNINGS
//
//ElementType array_to_type(const bp::numeric::array& a) {
//  return num_to_type(PyArray_DESCR(TP_ARRAY(a))->type_num);
//}
//
//size_t array_to_ndim(const bp::numeric::array& a) {
//  return PyArray_NDIM(TP_ARRAY(a));
//}
//
//#define TP_DESCR(x) ((PyArray_Descr*)x.ptr())
//
//dtype::dtype (bp::object dtype_like) {
//  PyArray_Descr* tmp = 0;
//  if (!PyArray_DescrConverter2(dtype_like.ptr(), &tmp)) {
//    std::string dtype_str = bp::extract<std::string>(bp::str(dtype_like));
//    PYTHON_ERROR(TypeError, "cannot convert input dtype-like object (%s) to proper dtype", dtype_str.c_str());
//  }
//  bp::handle<> hdl(bp::borrowed((PyObject*)tmp));
//  m_self = bp::object(hdl);
//}
//
//dtype::dtype (PyArray_Descr* descr) {
//  if (descr) {
//    bp::handle<> hdl((PyObject*)descr); //< raises if NULL
//    m_self = bp::object(hdl);
//  }
//}
//
//dtype::dtype(int typenum) {
//  PyArray_Descr* tmp = PyArray_DescrFromType(typenum);
//  bp::handle<> hdl(bp::borrowed((PyObject*)tmp));
//  m_self = bp::object(hdl);
//}
//
//dtype::dtype(ElementType eltype) {
//  if (eltype != t_unknown) {
//    PyArray_Descr* tmp = PyArray_DescrFromType(type_to_num(eltype));
//    bp::handle<> hdl(bp::borrowed((PyObject*)tmp));
//    m_self = bp::object(hdl);
//  }
//}
//
//dtype::dtype(const dtype& other): m_self(other.m_self)
//{
//}
//
//dtype::dtype() {
//}
//
//dtype::~dtype() { }
//
//dtype& dtype::operator= (const dtype& other) {
//  m_self = other.m_self;
//  return *this;
//}
//
//bool dtype::has_native_byteorder() const {
//  return TPY_ISNONE(m_self)? false : (PyArray_EquivByteorders(TP_DESCR(m_self)->byteorder, NPY_NATIVE) || TP_DESCR(m_self)->elsize == 1);
//}
//
//bool dtype::has_type(ElementType _eltype) const {
//  return eltype() == _eltype;
//}
//
//ElementType dtype::eltype() const {
//  return TPY_ISNONE(m_self)? t_unknown :
//    num_to_type(TP_DESCR(m_self)->type_num);
//}
//
//int dtype::type_num() const {
//  return TPY_ISNONE(m_self)? -1 : TP_DESCR(m_self)->type_num;
//}
//
//bp::str dtype::str() const {
//  return bp::str(m_self);
//}
//
//std::string dtype::cxx_str() const {
//  return bp::extract<std::string>(this->str());
//}
//
//PyArray_Descr* dtype::descr() {
//  return (PyArray_Descr*)m_self.ptr();
//}
//
///****************************************************************************
// * Free methods                                                             *
// ****************************************************************************/
//
//void typeinfo_ndarray_ (const bp::object& o, typeinfo& i) {
//  PyArrayObject* npy = TP_ARRAY(o);
//  npy_intp strides[NPY_MAXDIMS];
//  for (int k=0; k<PyArray_NDIM(npy); ++k) strides[k] = PyArray_STRIDES(npy)[k]/PyArray_DESCR(npy)->elsize;
//  i.set<npy_intp>(num_to_type(PyArray_DESCR(npy)->type_num), PyArray_NDIM(npy), PyArray_DIMS(npy), strides);
//}
//
//void typeinfo_ndarray (const bp::object& o, typeinfo& i) {
//  if (!PyArray_Check(o.ptr())) {
//    throw std::runtime_error("invalid input: cannot extract typeinfo object from anything else than ndarray");
//  }
//  typeinfo_ndarray_(o, i);
//}
//
///**
// * This method emulates the behavior of PyArray_GetArrayParamsFromObject from
// * NumPy >= 1.6 and is used when compiling and liking against older versions of
// * NumPy.
// */
//static int _GetArrayParamsFromObject(PyObject* op,
//    PyArray_Descr* requested_dtype,
//    npy_bool writeable,
//    PyArray_Descr** out_dtype,
//    int* out_ndim,
//    npy_intp* out_dims,
//    PyArrayObject** out_arr,
//    PyObject*) {
//
//  if (PyArray_Check(op)) { //it is already an array, easy
//
//    PyArrayObject* arr = reinterpret_cast<PyArrayObject*>(op);
//
//    if (requested_dtype && !PyArray_EquivTypes(PyArray_DESCR(arr), requested_dtype)) {
//
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//      if (PyArray_CanCastTypeTo(PyArray_DESCR(arr), requested_dtype, NPY_SAFE_CASTING)) {
//#else
//      if (PyArray_CanCastTo(PyArray_DESCR(arr), requested_dtype)) {
//#endif
//        (*out_arr) = 0;
//        (*out_dtype) = PyArray_DESCR(arr);
//        (*out_ndim) = PyArray_NDIM(arr);
//        for (int i=0; i<PyArray_NDIM(arr); ++i)
//          out_dims[i] = PyArray_DIM(arr,i);
//        // we need to cast the array, write-ability will not hold...
//        return writeable? 1 : 0;
//      }
//
//      else {
//        return 1;
//      }
//
//    }
//
//    //if you get to this point, the types are equivalent or there was no type
//    (*out_arr) = (PyArrayObject*)PyArray_FromArray(arr, 0, 0);
//    (*out_dtype) = 0;
//    (*out_ndim) = 0;
//    return writeable? (!PyArray_ISWRITEABLE(arr)) : 0;
//
//  }
//
//  else { //it is not an array -- try a brute-force conversion
//
//    TDEBUG1("[non-optimal] using NumPy version < 1.6 requires we convert input data for convertibility check - compile against NumPy >= 1.6 to improve performance");
//    bp::handle<> hdl(bp::allow_null(PyArray_FromAny(op, requested_dtype, 0, 0, 0, 0)));
//    bp::object array(hdl);
//
//    if (TPY_ISNONE(array)) return 1;
//
//    //if the conversion worked, you can now fill in the parameters
//    (*out_arr) = 0;
//    (*out_dtype) = PyArray_DESCR(TP_ARRAY(array));
//    (*out_ndim) = PyArray_NDIM(TP_ARRAY(array));
//    for (int i=0; i<PyArray_NDIM(TP_ARRAY(array)); ++i)
//      out_dims[i] = PyArray_DIM(TP_ARRAY(array),i);
//
//    //in this mode, the resulting object will never be write-able.
//    return writeable? 1 : 0;
//
//  }
//
//  return 0; //turn-off c compiler warnings...
//
//}
//
//convert_t convertible(bp::object array_like, typeinfo& info,
//    bool writeable, bool behaved) {
//
//  int ndim = 0;
//  npy_intp dims[NPY_MAXDIMS];
//  PyArrayObject* arr = 0;
//  PyArray_Descr* dtype = 0;
//
//  int not_convertible =
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//    PyArray_GetArrayParamsFromObject
//#else
//    _GetArrayParamsFromObject
//#endif
//    (array_like.ptr(), //input object pointer
//     0,                //requested dtype (if need to enforce)
//     writeable,        //writeable?
//     &dtype,           //dtype assessment - borrowed
//     &ndim,            //assessed number of dimensions
//     dims,             //assessed shape
//     &arr,             //if obj_ptr is ndarray, return it here
//     0)                //context?
//    ;
//
//  if (not_convertible) return IMPOSSIBLE;
//
//  convert_t retval = BYREFERENCE;
//
//  if (arr) { //the passed object is an array
//
//    //checks behavior.
//    if (behaved && !PyArray_ISCARRAY_RO(arr)) retval = WITHARRAYCOPY;
//
//    info.set<npy_intp>(num_to_type(PyArray_DESCR(arr)->type_num),
//        PyArray_NDIM(arr), PyArray_DIMS(arr));
//
//    Py_XDECREF(arr);
//  }
//
//  else { //the passed object is not an array
//    info.set<npy_intp>(num_to_type(dtype->type_num), ndim, dims);
//    retval = WITHCOPY;
//  }
//
//  return retval;
//}
//
//convert_t convertible_to (bp::object array_like,
//    const typeinfo& info, bool writeable, bool behaved) {
//
//  dtype req_dtype(info.dtype);
//
//  int ndim = 0;
//  npy_intp dims[NPY_MAXDIMS];
//  PyArrayObject* arr = 0;
//  PyArray_Descr* dtype = 0;
//
//  int not_convertible =
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//    PyArray_GetArrayParamsFromObject
//#else
//    _GetArrayParamsFromObject
//#endif
//    (array_like.ptr(),  //input object pointer
//     req_dtype.descr(), //requested dtype (if need to enforce)
//     writeable,         //writeable?
//     &dtype,            //dtype assessment
//     &ndim,             //assessed number of dimensions
//     dims,              //assessed shape
//     &arr,              //if obj_ptr is ndarray, return it here
//     0)                 //context?
//    ;
//
//  if (not_convertible) return IMPOSSIBLE;
//
//  convert_t retval = BYREFERENCE;
//
//  if (arr) { //the passed object is an array -- check compatibility
//
//    if (info.nd) { //check number of dimensions and shape, if needs to
//      if (PyArray_NDIM(arr) != (int)info.nd) {
//        Py_XDECREF(arr);
//        return IMPOSSIBLE;
//      }
//      if (info.has_valid_shape()) {
//        for (size_t i=0; i<info.nd; ++i)
//          if ((int)info.shape[i] != PyArray_DIM(arr,i)) {
//            Py_XDECREF(arr);
//            return IMPOSSIBLE;
//          }
//      }
//    }
//
//    //checks dtype
//    if (PyArray_DESCR(arr)->type_num != req_dtype.descr()->type_num) {
//      //well... we would need to cast, but maybe possible - check with NumPy
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//      if (PyArray_CanCastTypeTo(PyArray_DESCR(arr), req_dtype.descr(), NPY_SAFE_CASTING)) {
//#else
//      if (PyArray_CanCastTo(PyArray_DESCR(arr), req_dtype.descr())) {
//#endif
//        retval = WITHARRAYCOPY;
//      }
//      else {
//        //we cannot cast from current to desired type, sorry...
//        Py_XDECREF(arr);
//        return IMPOSSIBLE;
//      }
//    }
//
//    //checks behavior.
//    if (behaved) {
//      if (!PyArray_ISCARRAY_RO(arr)) retval = WITHARRAYCOPY;
//    }
//
//    Py_XDECREF(arr);
//
//  }
//
//  else { //the passed object is not an array
//
//    retval = WITHCOPY;
//
//    if (info.nd) { //check number of dimensions and shape
//      if (ndim != (int)info.nd) return IMPOSSIBLE;
//      for (size_t i=0; i<info.nd; ++i)
//        if (info.shape[i] &&
//            (int)info.shape[i] != dims[i]) return IMPOSSIBLE;
//    }
//
//    //checks dtype
//    if (dtype->type_num != req_dtype.descr()->type_num) {
//      //well... we would need to cast, but maybe possible - check with NumPy
//      if (
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//          !PyArray_CanCastTypeTo(dtype, req_dtype.descr(), NPY_SAFE_CASTING)
//#else
//          !PyArray_CanCastTo(dtype, req_dtype.descr())
//#endif
//         )
//      {
//        return IMPOSSIBLE;
//      }
//
//    }
//
//  }
//
//  return retval;
//}
//
//convert_t convertible_to(bp::object array_like, bp::object dtype_like,
//    bool writeable, bool behaved) {
//
//  dtype req_dtype(dtype_like);
//
//  int ndim = 0;
//  npy_intp dims[NPY_MAXDIMS];
//  PyArrayObject* arr = 0;
//  PyArray_Descr* dtype = 0;
//
//  int not_convertible =
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//    PyArray_GetArrayParamsFromObject
//#else
//    _GetArrayParamsFromObject
//#endif
//    (array_like.ptr(),  //input object pointer
//     req_dtype.descr(), //requested dtype (if need to enforce)
//     writeable,         //writeable?
//     &dtype,            //dtype assessment
//     &ndim,             //assessed number of dimensions
//     dims,              //assessed shape
//     &arr,              //if obj_ptr is ndarray, return it here
//     0)                 //context?
//    ;
//
//  if (not_convertible) return IMPOSSIBLE;
//
//  convert_t retval = BYREFERENCE;
//
//  if (arr) { //the passed object is an array -- check compatibility
//
//    //checks dtype
//    if (PyArray_DESCR(arr)->type_num != req_dtype.descr()->type_num) {
//      //well... we would need to cast, but maybe possible - check with NumPy
//
//      if (
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//          PyArray_CanCastTypeTo(PyArray_DESCR(arr), req_dtype.descr(), NPY_SAFE_CASTING)
//#else
//          PyArray_CanCastTo(PyArray_DESCR(arr), req_dtype.descr())
//#endif
//         )
//      {
//        retval = WITHARRAYCOPY;
//      }
//      else {
//        //we cannot cast from current to desired type, sorry...
//        Py_XDECREF(arr);
//        return IMPOSSIBLE;
//      }
//    }
//
//    //checks behavior.
//    if (behaved) {
//      if (!PyArray_ISCARRAY_RO(arr)) retval = WITHARRAYCOPY;
//    }
//
//    Py_XDECREF(arr);
//
//  }
//
//  else { //the passed object is not an array
//
//    retval = WITHCOPY;
//
//    //checks dtype
//    if (dtype->type_num != req_dtype.descr()->type_num) {
//      //well... we would need to cast, but maybe possible - check with NumPy
//      if (
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//          !PyArray_CanCastTypeTo(dtype, req_dtype.descr(), NPY_SAFE_CASTING)
//#else
//          !PyArray_CanCastTo(dtype, req_dtype.descr())
//#endif
//         )
//      {
//        return IMPOSSIBLE;
//      }
//
//    }
//
//  }
//
//  return retval;
//}
//
//convert_t convertible_to(bp::object array_like, bool writeable,
//    bool behaved) {
//
//  int ndim = 0;
//  npy_intp dims[NPY_MAXDIMS];
//  PyArrayObject* arr = 0;
//  PyArray_Descr* dtype = 0;
//
//  int not_convertible =
//#if NPY_FEATURE_VERSION >= NUMPY16_API /* NumPy C-API version >= 1.6 */
//    PyArray_GetArrayParamsFromObject
//#else
//    _GetArrayParamsFromObject
//#endif
//    (array_like.ptr(), //input object pointer
//     0,                //requested dtype (if need to enforce)
//     writeable,        //writeable?
//     &dtype,           //dtype assessment
//     &ndim,            //assessed number of dimensions
//     dims,             //assessed shape
//     &arr,             //if obj_ptr is ndarray, return it here
//     0)                //context?
//    ;
//
//  if (not_convertible) return IMPOSSIBLE;
//
//  convert_t retval = BYREFERENCE;
//
//  if (arr) { //the passed object is an array -- check compatibility
//
//    //checks behavior.
//    if (behaved) {
//      if (!PyArray_ISCARRAY_RO(arr)) retval = WITHARRAYCOPY;
//    }
//
//    Py_XDECREF(arr);
//
//  }
//
//  else { //the passed object is not an array
//
//     retval = WITHCOPY;
//
//  }
//
//  return retval;
//}
//
///***************************************************************************
// * Ndarray (PyArrayObject) manipulations                                   *
// ***************************************************************************/
//
///**
// * Returns either a reference or a copy of the given array_like object,
// * depending on the following requirements for referral:
// *
// * 0. The pointed object is a numpy.ndarray
// * 1. The array type description type_num matches
// * 2. The array is C-style, contiguous and aligned
// */
//static bp::object try_refer_ndarray (bp::object array_like,
//    bp::object dtype_like) {
//
//  PyArrayObject* candidate = TP_ARRAY(array_like);
//  PyArray_Descr* req_dtype = 0;
//  PyArray_DescrConverter2(dtype_like.ptr(), &req_dtype); //new ref!
//
//  bool can_refer = true; //< flags a copy of the data
//
//  if (!PyArray_Check((PyObject*)candidate)) can_refer = false;
//
//  if (can_refer && !PyArray_ISCARRAY_RO(candidate)) can_refer = false;
//
//  if (can_refer) {
//    PyObject* tmp = PyArray_FromArray(candidate, 0, 0);
//    bp::handle<> hdl(tmp); //< raises if NULL
//    bp::object retval(hdl);
//    return retval;
//  }
//
//  //copy
//  TDEBUG1("[non-optimal] copying array-like object - cannot refer");
//  PyObject* _ptr = (PyObject*)candidate;
//#if NPY_FEATURE_VERSION > NUMPY16_API /* NumPy C-API version > 1.6 */
//  int flags = NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ENSURECOPY|NPY_ARRAY_ENSUREARRAY;
//#else
//  int flags = NPY_C_CONTIGUOUS|NPY_ENSURECOPY|NPY_ENSUREARRAY;
//#endif
//  PyObject* tmp = PyArray_FromAny(_ptr, req_dtype, 0, 0, flags, 0);
//  bp::handle<> hdl(tmp); //< raises if NULL
//  bp::object retval(hdl);
//  return retval;
//
//}
//
//static void derefer_ndarray (PyArrayObject* array) {
//  Py_XDECREF(array);
//}
//
//static boost::shared_ptr<void> shared_from_ndarray (bp::object& o) {
//  boost::shared_ptr<PyArrayObject> cache(TP_ARRAY(o),
//      std::ptr_fun(derefer_ndarray));
//  Py_XINCREF(TP_OBJECT(o)); ///< makes sure it outlives this scope!
//  return cache; //casts to b::shared_ptr<void>
//}
//
//py_array::py_array(bp::object o, bp::object _dtype):
//  m_is_numpy(true)
//{
//  if (TPY_ISNONE(o)) PYTHON_ERROR(TypeError, "You cannot pass 'None' as input parameter to C++-bound bob methods that expect NumPy ndarrays (or blitz::Array<T,N>'s). Double-check your input!");
//  bp::object mine = try_refer_ndarray(o, _dtype);
//
//  //captures data from a numeric::array
//  typeinfo_ndarray_(mine, m_type);
//
//  //transforms the from bp ref counting to boost::shared_ptr<void>
//  m_data = shared_from_ndarray(mine);
//
//  //set-up the C-style pointer to this data
//  m_ptr = static_cast<void*>(PyArray_DATA(TP_ARRAY(mine)));
//}
//
//py_array::py_array(const interface& other) {
//  set(other);
//}
//
//py_array::py_array(boost::shared_ptr<interface> other) {
//  set(other);
//}
//
//py_array::py_array(const typeinfo& info) {
//  set(info);
//}
//
//py_array::~py_array() {
//}
//
///**
// * Wrap a C-style pointer with a PyArrayObject
// */
//static bp::object wrap_data (void* data, const typeinfo& ti,
//    bool writeable=true) {
//
//  npy_intp shape[NPY_MAXDIMS];
//  npy_intp stride[NPY_MAXDIMS];
//  for (size_t k=0; k<ti.nd; ++k) {
//    shape[k] = ti.shape[k];
//    stride[k] = ti.item_size()*ti.stride[k];
//  }
//  PyObject* tmp = PyArray_New(&PyArray_Type, ti.nd,
//        &shape[0], type_to_num(ti.dtype), &stride[0], data, 0,
//#if NPY_FEATURE_VERSION > NUMPY16_API /* NumPy C-API version > 1.6 */
//        writeable? NPY_ARRAY_CARRAY : NPY_ARRAY_CARRAY_RO
//#else
//        writeable? NPY_CARRAY : NPY_CARRAY_RO
//#endif
//        ,0);
//
//  bp::handle<> hdl(tmp);
//  bp::object retval(hdl);
//  return retval;
//}
//
///**
// * New wrapper of the array
// */
//static bp::object wrap_ndarray (const bp::object& a) {
//  PyObject* tmp = PyArray_FromArray(TP_ARRAY(a), 0, 0);
//  bp::handle<> hdl(tmp); //< raises if NULL
//  bp::object retval(hdl);
//  return retval;
//}
//
///**
// * Creates a new array from specifications
// */
//static bp::object make_ndarray(int nd, npy_intp* dims, int type) {
//  PyObject* tmp = PyArray_SimpleNew(nd, dims, type);
//  bp::handle<> hdl(tmp); //< raises if NULL
//  bp::object retval(hdl);
//  return retval;
//}
//
///**
// * New copy of the array from another array
// */
//static bp::object copy_array (const bp::object& array) {
//  PyArrayObject* _p = TP_ARRAY(array);
//  bp::object retval = make_ndarray(PyArray_NDIM(_p), PyArray_DIMS(_p), PyArray_DESCR(_p)->type_num);
//  PyArray_CopyInto(TP_ARRAY(retval), TP_ARRAY(array));
//  return retval;
//}
//
///**
// * Copies a data pointer and type into a new numpy array.
// */
//static bp::object copy_data (const void* data, const typeinfo& ti) {
//  bp::object wrapped = wrap_data(const_cast<void*>(data), ti);
//  bp::object retval = copy_array (wrapped);
//  return retval;
//}
//
//void py_array::set(const interface& other) {
//  TDEBUG1("[non-optimal] buffer copying operation being performed for "
//      << other.type().str());
//
//  //performs a copy of the data into a numpy array
//  bp::object mine = copy_data(other.ptr(), m_type);
//
//  //captures data from a numeric::array
//  typeinfo_ndarray_(mine, m_type);
//
//  //transforms the from bp ref counting to boost::shared_ptr<void>
//  m_data = shared_from_ndarray(mine);
//
//  //set-up the C-style pointer to this data
//  m_ptr = static_cast<void*>(PyArray_DATA(TP_ARRAY(mine)));
//
//  m_is_numpy = true;
//}
//
//void py_array::set(boost::shared_ptr<interface> other) {
//  m_type = other->type();
//  m_is_numpy = false;
//  m_ptr = other->ptr();
//  m_data = other->owner();
//}
//
///**
// * Creates a new numpy array from a bob::io::typeinfo object.
// */
//static bp::object new_from_type (const typeinfo& ti) {
//  npy_intp shape[NPY_MAXDIMS];
//  npy_intp stride[NPY_MAXDIMS];
//  for (size_t k=0; k<ti.nd; ++k) {
//    shape[k] = ti.shape[k];
//    stride[k] = ti.item_size()*ti.stride[k];
//  }
//  PyObject* tmp = PyArray_New(&PyArray_Type, ti.nd, &shape[0],
//      type_to_num(ti.dtype), &stride[0], 0, 0, 0, 0);
//  bp::handle<> hdl(tmp); //< raises if NULL
//  bp::object retval(hdl);
//  return retval;
//}
//
//void py_array::set (const typeinfo& req) {
//  if (m_type.is_compatible(req)) return; ///< nothing to do!
//
//  TDEBUG1("[non-optimal?] buffer re-size being performed from " << m_type.str()
//      << " to " << req.str());
//
//  bp::object mine = new_from_type(req);
//
//  //captures data from a numeric::array
//  typeinfo_ndarray_(mine, m_type);
//
//  //transforms the from bp ref counting to boost::shared_ptr<void>
//  m_data = shared_from_ndarray(mine);
//
//  //set-up the C-style pointer to this data
//  m_ptr = static_cast<void*>(PyArray_DATA(TP_ARRAY(mine)));
//
//  m_is_numpy = true;
//}
//
//bp::object py_array::copy(const bp::object& dtype) {
//  return copy_data(m_ptr, m_type);
//}
//
///**
// * Gets a read-only reference to a certain data. This recipe was originally
// * posted here:
// * http://blog.enthought.com/python/numpy-arrays-with-pre-allocated-memory/
// *
// * But a better allocation strategy (that actually works) is posted here:
// * http://stackoverflow.com/questions/2924827/numpy-array-c-api
// */
//// PyCapsule is only available on Python2.7 or Python3.1 and up
//#if ((PY_VERSION_HEX <  0x02070000) \
//    || ((PY_VERSION_HEX >= 0x03000000) \
//      && (PY_VERSION_HEX <  0x03010000)) )
//static void DeleteSharedPointer (void* ptr) {
//  typedef boost::shared_ptr<const void> type;
//  delete static_cast<type*>(ptr);
//}
//#else
//static void DeleteSharedPointer (PyObject* ptr) {
//  typedef boost::shared_ptr<const void> type;
//  delete static_cast<type*>(PyCapsule_GetPointer(ptr, NULL));
//}
//#endif
//
//static bp::object make_readonly (const void* data, const typeinfo& ti,
//    boost::shared_ptr<const void> owner) {
//
//  bp::object retval = wrap_data(const_cast<void*>(data), ti, false);
//
//  //creates the shared pointer deallocator
//  boost::shared_ptr<const void>* ptr = new boost::shared_ptr<const void>(owner);
//#if ((PY_VERSION_HEX <  0x02070000) \
//    || ((PY_VERSION_HEX >= 0x03000000) \
//      && (PY_VERSION_HEX <  0x03010000)) )
//  PyObject* py_sharedptr = PyCObject_FromVoidPtr(ptr, DeleteSharedPointer);
//#else
//  PyObject* py_sharedptr = PyCapsule_New((void*)ptr, NULL, DeleteSharedPointer);
//#endif
//
//  if (!py_sharedptr) {
//    PYTHON_ERROR(RuntimeError, "could not allocate space for deallocation object in read-only array::interface wrapping");
//  }
//
//# if NPY_FEATURE_VERSION >= NUMPY17_API /* NumPy C-API version >= 1.7 */
//  PyArray_SetBaseObject(TP_ARRAY(retval), py_sharedptr);
//# else
//  TP_ARRAY(retval)->base = py_sharedptr;
//# endif
//
//  return retval;
//}
//
//bp::object py_array::pyobject() {
//  if (m_is_numpy) {
//    bp::handle<> hdl(bp::borrowed(boost::static_pointer_cast<PyObject>(m_data).get()));
//    bp::object mine(hdl);
//    return mine;
//  }
//
//  //if you really want, I can wrap it for you, but in this case I'll make it
//  //read-only and will associate the object deletion to my own data pointer.
//  return make_readonly(m_ptr, m_type, m_data);
//}
//
//bool py_array::is_writeable() const {
//  return (!m_is_numpy || PyArray_ISWRITEABLE(boost::static_pointer_cast<PyArrayObject>(m_data).get()));
//}
//
//ndarray::ndarray(bp::object array_like, bp::object dtype_like)
//  : px(new py_array(array_like, dtype_like)) {
//}
//
//ndarray::ndarray(bp::object array_like)
//  : px(new py_array(array_like, bp::object())) {
//  }
//
//ndarray::ndarray(const typeinfo& info)
//  : px(new py_array(info)) {
//  }
//
//ndarray::~ndarray() { }
//
//const typeinfo& ndarray::type() const {
//  return px->type();
//}
//
//bp::object ndarray::self() { return px->pyobject(); }
//
//const_ndarray::const_ndarray(bp::object array_like)
//  : ndarray(array_like) {
//  }
//
//const_ndarray::~const_ndarray() { }
//
//}}
//
