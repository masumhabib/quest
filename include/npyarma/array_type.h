/**
 * @file bob/core/array_type.h
 * @date Sat Apr 9 18:10:10 2011 +0200
 * @author Laurent El Shafey <Laurent.El-Shafey@idiap.ch>
 *
 * @brief This file contains information about the supported arrays
 *
 * Copyright (C) 2011-2013 Idiap Research Institute, Martigny, Switzerland
 */

#ifndef BOB_CORE_ARRAY_TYPE_H
#define BOB_CORE_ARRAY_TYPE_H

#include <stdint.h>
#include <cstdlib>
#include <complex>

/**
 * @ingroup CORE_ARRAY
 * @brief This macro defines the maximum number of dimensions supported by bob.
 * A variable in the bob.core.array namespace is created from this macro
 * receiving the same value. Use that variable on your programs, or this macro
 * on your preprocessor code.
 */
#define BOB_MAX_DIM 4

namespace bob {
  namespace core { namespace array {
    /**
     * @ingroup CORE_ARRAY
     * @{
     */

    /**
     * @brief Enumeration of the supported type for multidimensional arrays
     * @warning float128 and complex256 are defined but currently not
     * supported
     */
    typedef enum ElementType {
      t_unknown=0,
      t_bool=1,
      t_int8=2,
      t_int16=3,
      t_int32=4,
      t_int64=5,
      t_uint8=6,
      t_uint16=7,
      t_uint32=8,
      t_uint64=9,
      t_float32=10,
      t_float64=11,
      t_float128=12,
      t_complex64=13,
      t_complex128=14,
      t_complex256=15
    } ElementType;

    /**
     * @brief Maximum number of supported dimensions for multidimensional
     * arrays.
     */
    const size_t N_MAX_DIMENSIONS_ARRAY = BOB_MAX_DIM;

    /**
     * @brief These are some type to element type conversions
     */
    template<typename T> ElementType getElementType() {
      return t_unknown;
    }

    /**
     * @brief Some specializations that convert type to element type.
     */
    template<> inline ElementType getElementType<bool>() { return t_bool; }
    template<> inline ElementType getElementType<int8_t>() { return t_int8; }
    template<> inline ElementType getElementType<int16_t>()
    { return t_int16; }
    template<> inline ElementType getElementType<int32_t>()
    { return t_int32; }
    template<> inline ElementType getElementType<int64_t>()
    { return t_int64; }
    template<> inline ElementType getElementType<uint8_t>()
    { return t_uint8; }
    template<> inline ElementType getElementType<uint16_t>()
    { return t_uint16; }
    template<> inline ElementType getElementType<uint32_t>()
    { return t_uint32; }
    template<> inline ElementType getElementType<uint64_t>()
    { return t_uint64; }
    template<> inline ElementType getElementType<float>()
    { return t_float32; }
    template<> inline ElementType getElementType<double>()
    { return t_float64; }
    template<> inline ElementType getElementType<long double>()
    { return t_float128; }
    template<> inline ElementType getElementType<std::complex<float> >()
    { return t_complex64; }
    template<> inline ElementType getElementType<std::complex<double> >()
    { return t_complex128; }
    template<> inline ElementType getElementType<std::complex<long double> >()
    { return t_complex256; }

    /**
     * @brief These are some type to element size conversions
     */
    template<typename T> size_t getElementSize() {
      return 0;
    }

    /**
     * @brief Some specializations that convert the types we handle properly
     */
    template<> inline size_t getElementSize<bool>() { return sizeof(bool); }
    template<> inline size_t getElementSize<int8_t>()
    { return sizeof(int8_t); }
    template<> inline size_t getElementSize<int16_t>()
    { return sizeof(int16_t); }
    template<> inline size_t getElementSize<int32_t>()
    { return sizeof(int32_t); }
    template<> inline size_t getElementSize<int64_t>()
    { return sizeof(int64_t); }
    template<> inline size_t getElementSize<uint8_t>()
    { return sizeof(uint8_t); }
    template<> inline size_t getElementSize<uint16_t>()
    { return sizeof(uint16_t); }
    template<> inline size_t getElementSize<uint32_t>()
    { return sizeof(uint32_t); }
    template<> inline size_t getElementSize<uint64_t>()
    { return sizeof(uint64_t); }
    template<> inline size_t getElementSize<float>()
    { return sizeof(float); }
    template<> inline size_t getElementSize<double>()
    { return sizeof(double); }
    template<> inline size_t getElementSize<long double>()
    { return sizeof(long double); }
    template<> inline size_t getElementSize<std::complex<float> >()
    { return sizeof(std::complex<float>); }
    template<> inline size_t getElementSize<std::complex<double> >()
    { return sizeof(std::complex<double>); }
    template<> inline size_t getElementSize<std::complex<long double> >()
    { return sizeof(std::complex<long double>); }

    /**
     * @brief Returns the type size given the enumeration
     */
    size_t getElementSize(ElementType t);

    /**
     * @brief Gets a string representation of an element type value
     */
    const char* stringize(ElementType t);

    /**
     * @brief Equivalent to call stringize() on the result of 
     * getElementType<T>().
     */
    template<typename T> const char* stringize() {
      return stringize(getElementType<T>());
    }

    /**
     * @brief Returns the ElementType given the string representation
     */
    ElementType unstringize(const char* name);

    /**
     * @}
     */
  }}
}

#endif /* BOB_CORE_ARRAY_TYPE_H */
