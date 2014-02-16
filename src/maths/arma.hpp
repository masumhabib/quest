/* 
 * File:   arma.hpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 8, 2014, 8:12 PM
 */

#ifndef ARMA_HPP
#define	ARMA_HPP

#include <complex>
#include <armadillo>

namespace maths{
namespace armadillo{

typedef std::complex<double>      dcmplx;   //!< Double precision complex
typedef std::complex<float>       fcmplx;   //!< Single precision complex

using arma::Mat;
using arma::Col;
using arma::Row;
using arma::span;

using arma::field;
using arma::endr;
using arma::inv;
using arma::zeros;
using arma::ones;
using arma::linspace;
using arma::eye;
using arma::is_finite;
using arma::wall_clock;
namespace fill = arma::fill;


typedef arma::Col<double>         vec;       //!< Double precision column vector.
typedef arma::Col<dcmplx>         cxvec;     //!< Double precision complex column vector.
typedef arma::Col<float>          fvec;      //!< Double precision column vector.
typedef arma::Col<fcmplx>         cxfvec;    //!< Double precision complex column vector.
typedef arma::Col<int>            ivec;      //!< Signed integer column vector.
typedef arma::Col<uint>           uvec;      //!< Unsigned integer column vector.

typedef arma::Col<double>         col;       //!< Double precision column vector.
typedef arma::Col<dcmplx>         cxcol;     //!< Double precision complex column vector.
typedef arma::Col<float>          fcol;      //!< Double precision column vector.
typedef arma::Col<fcmplx>         cxfcol;    //!< Double precision complex column vector.
typedef arma::Col<int>            icol;      //!< Signed integer column vector.
typedef arma::Col<uint>           ucol;      //!< Unsigned integer column vector.

typedef arma::Row<double>         row;       //!< Double precision row vector.
typedef arma::Row<dcmplx>         cxrow;     //!< Double precision complex row vector.
typedef arma::Row<float>          frow;      //!< Double precision row vector.
typedef arma::Row<fcmplx>         cxfrow;    //!< Double precision complex row vector.
typedef arma::Row<int>            irow;      //!< Signed integer column vector.
typedef arma::Row<uint>           urow;      //!< Unsigned integer column vector.

typedef arma::Mat<double>         mat;       //!< Double precision matrix.
typedef arma::Mat<dcmplx>         cxmat;     //!< Double precision complex matrix.
typedef arma::Mat<float>          fmat;      //!< Double precision matrix.
typedef arma::Mat<fcmplx>         cxfmat;    //!< Double precision complex matrix.
typedef arma::Mat<int>            imat;      //!< Signed integer matrix.
typedef arma::Mat<uint>           umat;      //!< Unsigned integer matrix.

typedef arma::Cube<double>         cube;       //!< Double precision matrix.
typedef arma::Cube<dcmplx>         cxcube;     //!< Double precision complex matrix.
typedef arma::Cube<float>          fcube;      //!< Double precision matrix.
typedef arma::Cube<fcmplx>         cxfcube;    //!< Double precision complex matrix.
typedef arma::Cube<int>            icube;      //!< Signed integer matrix.
typedef arma::Cube<uint>           ucube;      //!< Unsigned integer matrix.

typedef arma::Mat<double>::fixed<2,2>mat22;     //!< 2x2 Double precision matrix.
typedef cxmat::fixed<2,2>          cxmat22;   //!< 2x2 Double precision complex matrix.
typedef fmat::fixed<2,2>           fmat22;    //!< 2x2 Double precision matrix.
typedef cxfmat::fixed<2,2>         cxfmat22;  //!< 2x2 Double precision complex matrix.
typedef imat::fixed<2,2>           imat22;    //!< 2x2 Signed integer matrix.
typedef umat::fixed<2,2>           umat22;    //!< 2x2 Unsigned integer matrix.

typedef row::fixed<2>              row2;      //!< 2 element double precision row vector.
typedef row::fixed<3>              row3;      //!< 3 element double precision row vector.

}
}


#endif	/* ARMA_HPP */

