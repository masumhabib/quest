/* 
 * File:   ConsoleProgressBar.h
 * Author: masum
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
    static const float        mPercentPerDot    = 100.0/mnTotalDots;
    
    bool            mLivePercent;
    unsigned int    mnDots;
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
    
    virtual ~ConsoleProgressBar();
    
    ConsoleProgressBar& operator+= (unsigned long newCount);
    friend ConsoleProgressBar operator+ (ConsoleProgressBar lhs, unsigned long newCount);
    ConsoleProgressBar& operator++();
    ConsoleProgressBar operator++(int);
    
    void complete();
    void start();
    
       
protected:
    inline unsigned long convertCountToDots();
    virtual void draw();
    virtual void step(unsigned long count);
private:

};

}
#endif	/* CONSOLEPROGRESSBAR_H */

