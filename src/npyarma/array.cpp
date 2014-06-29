/**
 * @file core/cxx/array.cc
 * @date Tue Nov 8 15:34:31 2011 +0100
 * @author Andre Anjos <andre.anjos@idiap.ch>
 *
 * @brief Some buffer stuff
 *
 * Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland
 */

//#include <boost/format.hpp>
//#include "array.h"
//
//namespace qmicad{ namespace python{
//typeinfo::typeinfo():
//  dtype(t_unknown),
//  nd(0)
//{
//}
//
//typeinfo::typeinfo(const typeinfo& other): 
//  dtype(other.dtype)
//{
//  set_shape(other.nd, other.shape);
//}
//
//typeinfo& typeinfo::operator= (const typeinfo& other) {
//  dtype = other.dtype;
//  set_shape(other.nd, other.shape);
//  return *this;
//}
//
//void typeinfo::reset() {
//  dtype = t_unknown;
//  nd = 0;
//}
//
//bool typeinfo::is_valid() const {
//  return (dtype != t_unknown) && (nd > 0) && (nd <= (BOB_MAX_DIM+1)) && has_valid_shape();
//}
//
//void typeinfo::update_strides() {
//  switch (nd) {
//    case 0:
//      return;
//    case 1:
//      stride[0] = 1;
//      return;
//    case 2:
//      stride[1] = 1;
//      stride[0] = shape[1];
//      return;
//    case 3:
//      stride[2] = 1;
//      stride[1] = shape[2];
//      stride[0] = shape[1]*shape[2];
//      return;
//    case 4:
//      stride[3] = 1;
//      stride[2] = shape[3];
//      stride[1] = shape[2]*shape[3];
//      stride[0] = shape[1]*shape[2]*shape[3];
//      return;
//    case 5:
//      stride[4] = 1;
//      stride[3] = shape[4];
//      stride[2] = shape[3]*shape[4];
//      stride[1] = shape[2]*shape[3]*shape[4];
//      stride[0] = shape[1]*shape[2]*shape[3]*shape[4];
//      return;
//    default:
//      break;
//  }
//  throw std::runtime_error("unsupported number of dimensions");
//}
//
//size_t typeinfo::size() const {
//  size_t retval = 1;
//  for (size_t k=0; k<nd; ++k) retval *= shape[k];
//  return retval;
//}
//
//size_t typeinfo::buffer_size() const {
//  return size()*getElementSize(dtype);
//}
//
//static bool same_shape(size_t nd, const size_t* s1, const size_t* s2) {
//  for (size_t k=0; k<nd; ++k) if (s1[k] != s2[k]) return false;
//  return true;
//}
//
//bool typeinfo::is_compatible(const typeinfo& other) const {
//  return (dtype == other.dtype) && (nd == other.nd) && same_shape(nd, shape, other.shape);
//}
//
//std::string typeinfo::str() const {
//  boost::format s("dtype: %s (%d); shape: [%s]; size: %d bytes");
//  size_t sz = 0;
//  size_t buf_sz = 0;
//  if (dtype != t_unknown) {
//    //otherwise it throws
//    sz = item_size();
//    buf_sz = buffer_size();
//  }
//  s % item_str() % sz; 
//  switch (nd) {
//    case 0:
//      s % "";
//      break;
//    case 1:
//      s % (boost::format("%d") % shape[0]).str();
//      break;
//    case 2:
//      s % (boost::format("%d,%d") % shape[0] % shape[1]).str();
//      break;
//    case 3:
//      s % (boost::format("%d,%d,%d") % shape[0] % shape[1] % shape[2]).str();
//      break;
//    case 4:
//      s % (boost::format("%d,%d,%d,%d") % shape[0] % shape[1] % shape[2] % shape[3]).str();
//      break;
//    default:
//      s % ">4 dimensions?";
//      break;
//  }
//  s % buf_sz;
//  return s.str();
//}
//
//void typeinfo::reset_shape() {
//  shape[0] = 0;
//}
//
//bool typeinfo::has_valid_shape() const {
//  return shape[0] != 0;
//}
//
//
//}}