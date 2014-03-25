/* 
 * File:   Printable.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 18, 2013, 4:27 PM
 * 
 * Description: Printable class. All class that want to print itself
 * extends this class and overloads the toString() method.
 * 
 */

#ifndef PRINTABLE_H
#define	PRINTABLE_H

#include "std.hpp"
#include "../utils/serialize.hpp"

/**
 * Description: Printable class. All class that want to print itself
 * extends this class and overloads the toString() method.
 */

namespace utils{
using namespace stds;

class Printable{
protected:
    string mTitle;          //!< Object title.
    string mPrefix;         //!< Object prefix.

public:
        
    Printable(const string &prefix = ""):mPrefix(prefix) {
    };
    virtual ~Printable(){};

    virtual void print() const{
        cout << *this;
    }
    
    virtual string toString() const { 
        stringstream out;
        if (mTitle.length()){
                out << mPrefix << mTitle;
            }        
        return out.str(); 
    };

    virtual bool   isEmpty()const { return false; };

    /* Dump the data to the stream */
    friend ostream& operator << (ostream & out, const Printable &p){
        if(!p.isEmpty()){
            out << p.toString();
        }
        return out;
    };
    
    
    void Prefix(const string& prefix){
        mPrefix = prefix;  
    };

    string Prefix(){
        return mPrefix;
    };

    void Title(const string &title){
        mTitle = title;
    }

    string Title(){
        return mTitle;
    }
    
//!< Serialization
protected:    
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & mTitle;
        ar & mPrefix;
    }

};
}
#endif	/* PRINTABLE_H */

