/* 
 * File:   cache.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 2:50 PM
 */

#ifndef CACHE_HPP
#define	CACHE_HPP

#include <armadillo>
#include <stdexcept>
#include "maths/constants.h"
#include "utils/myenums.hpp"

namespace utils{
using std::out_of_range;
using namespace maths::armadillo;


/*
 * A cache of elements stored in a vector.
 */
template <class T>
class Cache{
public:
    Cache(int begin, int end, bool cacheEnabled = true){  
        mCacheEnabled = cacheEnabled;
        // Forward moving or backward moving
        if (begin <= end){
            mBegin = begin;
            mEnd = end;

        }else{
            mBegin = end;
            mEnd = begin;
        }

        // Initally we do cache is empty
        mIt = mBegin - 1;

        // Size of the cache
        mLength = mEnd - mBegin + 1;

        if (cacheEnabled == true){
            mM.set_size(mLength);
        }else{
            mM.set_size(1);
        }
    };

    // resets the cache.
    void reset(){
        // reset iterator
        mIt = mBegin - 1;
        // reset members.
        if (mCacheEnabled == true){
            mM.set_size(mLength);
        }else{
            mM.set_size(1);
        }        
    }

    // () operator is the read only access.
    virtual const T& operator()(int it){
        return getAt(it);
    };

    // convert easy index to array index (starting from 0)
    int toArrayIndx(int it){
        return it - mBegin;
    };

    void enableCache(bool enable = false){
        if (mCacheEnabled != enable){        
            mCacheEnabled = enable;
            reset();
        }
    }
    
    int begin(){ return mBegin; };
    int end(){ return mEnd; };
    int length(){ return mLength; };
    int current(){ return mIt; };

protected:    
    T& getAt(int it){
        // within the range
        if (it <= mEnd && it >= mBegin){             
            if (mCacheEnabled == true){
                int ii = toArrayIndx(it);
                return mM(ii);
            }else{
                return mM(0);
            }
        }else{
            // out of range, return empty object
            throw out_of_range("In Cache::getAt(i), i out of range.");
        }
    };

public:
protected:
    field<T>    mM;
    int         mIt;
    int         mBegin;
    int         mEnd;
    int         mLength;

    bool      mCacheEnabled;
private:
    Cache();
};

/*
 * Specialization of cache for armadillo matrix
 */
template<typename T>
class MatCache:public Cache<Mat<T> > {
public:
    MatCache( int begin, int end, bool cacheEnabled = true):
        Cache<Mat<T> >(begin, end, cacheEnabled){
    };
    
    
protected:
    bool isStored(int it){
        bool result;
        // within the range
        if (it <= this->mEnd && it >= this->mBegin){             
            if (this->mCacheEnabled == true){
                int ii = this->toArrayIndx(it);
                if (this->mM(ii).empty()){
                    result = false;
                }else{
                    result = true;
                }
            }else{
                if (it == this->mIt){
                    result = true;
                }else{
                    result = false;
                }
            }
        }else{
            result = false;
        }

        return result;
    };  
    
private:
    MatCache();
};

typedef MatCache<dcmplx> CxMatCache;

}

#endif	/* CACHE_HPP */

