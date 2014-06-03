/* 
 * File:   ConsoleProgressBar.h
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 8, 2013, 11:35 AM
 */

#ifndef CONSOLEPROGRESSBAR_H
#define	CONSOLEPROGRESSBAR_H

#include "vout.h"

#include <iostream>
#include <string>
#include <sstream>

namespace utils{
using namespace stds;

class ConsoleProgressBar {
protected:
    static const unsigned int mnTotalDots       = 60;
    
    unsigned long   mCount;
    unsigned long   mExpectedCount;
    
    string          mPrefix;
    string          mSuffix;
    char            mSpace;
    char            mDot;
    char            mEnds;
    
    
public:
    ConsoleProgressBar(string prefix = "", unsigned long expectedCount = 100,
                       unsigned long count = 0, bool livePercent = false);
    
    ConsoleProgressBar& operator+= (unsigned long newCount);
    friend ConsoleProgressBar operator+ (ConsoleProgressBar lhs, unsigned long newCount);
    ConsoleProgressBar& operator++();
    ConsoleProgressBar operator++(int);
    
    void complete();
    void start();
    
       
protected:
    virtual void step(unsigned long count);
private:

};

}
#endif	/* CONSOLEPROGRESSBAR_H */

