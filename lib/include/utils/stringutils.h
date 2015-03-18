/* 
 * File:   stringutils.h
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 5, 2013, 4:11 PM
 */

#ifndef STRINGUTILS_H
#define	STRINGUTILS_H
    
#include <time.h>
#include <boost/lexical_cast.hpp>
#include "utils/std.hpp"

namespace utils{namespace strings{
using stds::string;
using stds::stringstream;
using boost::lexical_cast;

string trim(const string& s, const char* t = " \t\n\r\f\v");

inline double stod(const string& ss){
    return lexical_cast<double>(trim(ss));
}


inline int stoi(const string& ss){
    return lexical_cast<int>(trim(ss));
}

inline long stou(const string& ss){
    return lexical_cast<long>(trim(ss));
}  

inline string dtos(double x){
    return lexical_cast<string>(x);
}

inline string itos(int x){
    return lexical_cast<string>(x);
}

inline string ttos(double sec){
    long tm = long(sec);
    long days = tm/(3600*24);
    tm = tm % (3600*24);

    long hr = tm/3600;
    tm = tm % 3600;

    long min = tm/60;
    double second = long(tm%60) + (sec - long(sec));

    stringstream times;
    times.precision(2);
    if (days > 0){
        times << days << " day(s) ";
    }

    if (hr > 0){
        times << hr << " hour(s) ";
    }

    if (min > 0){
        times << min << " minute(s) ";
    }

    times << second << " second(s)";

    return times.str();
}

}}
#endif	/* STRINGUTILS_H */

