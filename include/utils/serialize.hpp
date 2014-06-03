/* 
 * File:   serialize.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 11, 2013, 5:34 PM
 * 
 * Serialization for armadillo classes. The serialization is used by
 * the Boost::MPI to transmit and receive objects over the network.
 * 
 */

#ifndef SERIALIZE_HPP
#define	SERIALIZE_HPP

#include "maths/svec.h"

#include <string>
#include <sstream>
#include <armadillo>
#include <boost/serialization/access.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/split_free.hpp>

#include <boost/serialization/complex.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/list.hpp>
#include <boost/serialization/map.hpp>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>



namespace boost {
namespace serialization {
using std::stringstream;
using std::string;
using std::pair;
using std::complex;
using arma::Mat;
using arma::Col;
using arma::Row;
using namespace maths::spvec;

/*
 * Serialization for stringstream.
 * =============================================================================
 */

template<class Archive>
void save(Archive& ar, const stringstream& ss, const unsigned int version){
    string s = ss.str();
    ar & s;
};

template<class Archive>
void load(Archive& ar, stringstream& ss, const unsigned int version){
    string s;
    ar & s;
    ss.str(s);
};

template<class Archive>
inline void serialize(Archive& ar, stringstream& ss, const unsigned int file_version){
    split_free(ar, ss, file_version); 
};


/*
 * Serialization for svec.
 * =============================================================================
 */

template<class Archive>
void serialize(Archive& ar, svec& sv, const unsigned int version){
    ar & sv(coord::X);
    ar & sv(coord::Y);
    ar & sv(coord::Z);
};

/*
 * Serialization for arma::Mat<T>: Armadillo matrices.
 * V2: Using raw data.
 * =============================================================================
 */
template<class Archive, class T>
void save(Archive& ar, const Mat<T>& mat, const unsigned int version){
    uint n_rows = mat.n_rows;
    uint n_cols = mat.n_cols;
    
     ar << n_rows;
     ar << n_cols;
    for(uint ir = 0; ir < n_rows; ++ ir){
        for(uint ic = 0; ic < n_cols; ++ ic){
            ar << mat(ir, ic);
        }
    }
};

template<class Archive, class T>
void load(Archive& ar, Mat<T>& mat, const unsigned int version){
    uint n_rows;
    uint n_cols;
    
    ar >> n_rows;
    ar >> n_cols;
    
    if (n_rows != mat.n_rows || n_cols != mat.n_cols){
        mat.set_size(n_rows, n_cols);
    }
     
    for(int ir = 0; ir < n_rows; ++ ir){
        for(int ic = 0; ic < n_cols; ++ ic){
            ar >> mat(ir, ic);
        }
    }
};

template<class Archive, class T>
inline void serialize(Archive& ar, Mat<T>& m, const unsigned int file_version){
    split_free(ar, m, file_version); 
};

/*
 * Serialization for arma::Col<T>: Armadillo column vector.
 * V2: using raw data.
 * =============================================================================
 */
template<class Archive, class T>
void save(Archive& ar, const Col<T>& col, const unsigned int version){
    int n_rows = col.n_rows;
    ar << n_rows;
    for(int it = 0; it < col.n_rows; ++ it){
        ar << col(it);
    }
};

template<class Archive, class T>
void load(Archive& ar, Col<T>& col, const unsigned int version){
    uint n_rows;
    ar >> n_rows;
    if (n_rows != col.n_rows){
        col.set_size(n_rows);
    }
    for(uint it = 0; it < n_rows; ++ it){
        ar >> col(it);
    }
};

template<class Archive, class T>
inline void serialize(Archive& ar, Col<T>& col, const unsigned int file_version){
    split_free(ar, col, file_version); 
};

/*
 * Serialization for arma::Row<T>: Armadillo row vector.
 * V2: using raw data.
 * =============================================================================
 */
template<class Archive, class T>
void save(Archive& ar, const Row<T>& row, const unsigned int version){
    uint n_cols = row.n_cols;
    ar << n_cols;
    for(uint it = 0; it < n_cols; ++ it){
        ar << row(it);
    }
};

template<class Archive, class T>
void load(Archive& ar, Row<T>& row, const unsigned int version){
    uint n_cols;
    ar >> n_cols;
    if (n_cols != row.n_cols){
        row.set_size(n_cols);
    }
    for(uint it = 0; it < n_cols; ++ it){
        ar >> row(it);
    }
};

template<class Archive, class T>
inline void serialize(Archive& ar, Row<T>& row, const unsigned int file_version){
    split_free(ar, row, file_version); 
};

}} 


#endif	/* SERIALIZE_HPP */

