/* 
 * File:   property.hpp
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on May 29, 2014, 11:09 AM
 */

#ifndef PROPERTY_HPP
#define	PROPERTY_HPP

/**
 * C# like property in C++. Copied from Wikipedia.
 */
template <typename T> class property {
    T value;
public:
    T & operator = (const T &i) {
        return value = i;
    }
    template <typename T2> T2 & operator = (const T2 &i) {
        return value = i;

        // This template class member function template serves the purpose to make
        // typing more strict. Assignment to this is only possible with exact identical
        // types.

        //T2 &guard = value;
        //throw guard; // Never reached.
    }
    operator T const & () const {
        return value;
    }
};

#endif	/* PROPERTY_HPP */

