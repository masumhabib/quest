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

#ifndef QBASE_H
#define	QBASE_H

#include <string>
#include <iostream>

/**
 * Base class of all QMICAD classes. All class that want to print itself
 * extends this class and overloads the toString() method.
 */

namespace qmicad{
using std::ostream;
using std::endl;
using std::string;

class Qbase{
protected:
    string mTitle;
    string mPrefix;

    Qbase(const string &prefix = ""):mPrefix(prefix) {
        mTitle = "Qbase";
    };
    virtual ~Qbase(){};

public:
    virtual string toString() const { return mTitle; };
    virtual bool   isEmpty()const { return false; };

    /* Dump the data to the stream */
    friend ostream& operator << (ostream & out, const Qbase &p){
        if(!p.isEmpty()){
            if (p.mTitle.length()){
                out << p.mPrefix << p.mTitle << ":" << endl;
            }
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
};
}
#endif	/* QBASE_H */

