/* 
 * File:   ConsoleProgressBar.h
 * Author: masum
 *
 * Created on April 8, 2013, 11:35 AM
 */

#ifndef CONSOLEPROGRESSBAR_H
#define	CONSOLEPROGRESSBAR_H

#include <iostream>
#include <string>
#include <sstream>

using namespace std;

class ConsoleProgressBar {
protected:
    static const unsigned int mnTotalDots = 50;
    static const unsigned int mPercentPerDot = 2;
    
    unsigned int mnDots;
    unsigned long mCount;
    unsigned long mExpectedCount;
    
    string mPrefix;
    string mSuffix;
    char mSpace;
    char mDot;
    char mEnds;
    
    
public:
    ConsoleProgressBar(unsigned long expectedCount = 100, 
                       unsigned long count = 0, 
                       string prefix = "", 
                       string suffix = "%");
    virtual ~ConsoleProgressBar();
    
    ConsoleProgressBar& operator+= (unsigned long newCount);
    friend ConsoleProgressBar operator+ (ConsoleProgressBar lhs, unsigned long newCount);
    ConsoleProgressBar& operator++();
    ConsoleProgressBar operator++(int);
    
    void complete();
    void start();
    
       
protected:
    inline unsigned long convertCountToDots();
    void draw();
    void step(unsigned long count);
private:

};

#endif	/* CONSOLEPROGRESSBAR_H */

