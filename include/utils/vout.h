/* 
 * File:   vout.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 16, 2014, 12:31 PM
 */

#ifndef VOUT_H
#define	VOUT_H

#include "std.hpp"

namespace utils{
namespace stds{


/**
 * The verbosity structure handles the verbosity of the program.
 */

struct verbosity{
    verbosity(int verb):verb(verb){
    }
    
    verbosity(){
    }
    
    friend bool operator == (const verbosity &lhs, const verbosity &rhs){
        return lhs.verb == rhs.verb;
    }
    
    friend bool operator >= (const verbosity &lhs, const verbosity &rhs){
        return lhs.verb >= rhs.verb;
    }

    friend bool operator <= (const verbosity &lhs, const verbosity &rhs){
        return !operator <= (lhs, rhs);
    }

    int verb;

};

/**
 * Static verbosity levels.
 */
extern const verbosity  vquiet;
extern const verbosity  vnormal;
extern const verbosity  vdebug;
extern const verbosity  vdump;


/*
 * Verbosity controlled stream.
 */
class vostream{
public:
    vostream(ostream& out, const verbosity &verb = vnormal):
        mout(out), mAppVerb(verb), mCurrVerb(verb)
    {
        mMyId = 0;
        mPrintersId = 0;
    };
    
    template<typename T>
    vostream& operator<<(const T& v) {
        if (mAppVerb >= mCurrVerb && mPrintersId == mMyId){
            mout << v;
        }
        return *this;
    }
    
    vostream& operator<<(const verbosity& v) {
        mCurrVerb = v;
        return *this;
    }

    // For endl and such to work.
    ostream& operator<< ( ostream& (*pf)(ostream&)) {
        if (mAppVerb >= mCurrVerb && mPrintersId == mMyId){
            return pf(this->mout);
        }
        
        return this->mout;
    };
    
    void myId(int id){ mMyId = id; };
    int  myId() const { return mMyId; };
    void printersId(int id){ mPrintersId = id; };
    void appVerbosity(const verbosity &verb){ mAppVerb = verb; };
    
protected:
    ostream                   &mout;
    verbosity                 mAppVerb;         //!< Verbosity level of this program.
    verbosity                 mCurrVerb;        //!< Current verbosity level.
    int                       mMyId;            //!< ID of this process.
    int                       mPrintersId;      //!< ID of the printer process.
};

extern vostream         dout;
extern vostream         vout;

extern const string     dbg;

}
}
#endif	/* VOUT_H */

