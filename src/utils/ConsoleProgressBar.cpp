/* 
 * File:   ConsoleProgressBar.cpp
 * Author: masum
 * 
 * Created on April 8, 2013, 11:35 AM
 */

#include <ios>

#include "ConsoleProgressBar.h"    

ConsoleProgressBar::ConsoleProgressBar(unsigned long expectedCount, 
    unsigned long count, string prefix, string suffix):
mExpectedCount(expectedCount),
mCount(count < expectedCount ? count:expectedCount),
mPrefix(prefix),
mSuffix(suffix){
    
    mSpace = ' ';
    mDot   = '=';
    mEnds  = '|';

    mnDots = convertCountToDots();
}

ConsoleProgressBar::~ConsoleProgressBar() {
}

void ConsoleProgressBar::draw(){
    string bar;
    stringstream ssbar;
    
    bar.insert(bar.length(), mnDots, mDot);
    bar.insert(bar.length(), mnTotalDots - mnDots, mSpace);

    ssbar << mPrefix << "[";
    ssbar.width(3);
    ssbar << mnDots*mPercentPerDot;
    ssbar << mSuffix;
    ssbar << "] " << mEnds;
    ssbar << bar << mEnds;
    
    cout << ssbar.str() << endl;
    cout << "\033[F";
}

void ConsoleProgressBar::step(unsigned long count){
    
    mCount += count;
    if(mCount > mExpectedCount){
        mCount = mExpectedCount;
    }
    
    unsigned int newDots = convertCountToDots();
    if (newDots > mnDots){
        mnDots = newDots;
        draw();
    }    
}

unsigned long ConsoleProgressBar::convertCountToDots(){
    return (mCount*mnTotalDots)/mExpectedCount;
}

ConsoleProgressBar& ConsoleProgressBar::operator+= (unsigned long newCount){

    step(newCount);
    return *this;
}

ConsoleProgressBar operator+ (ConsoleProgressBar lhs, 
        unsigned long newCount){
    lhs += newCount;
    return lhs;
}

ConsoleProgressBar& ConsoleProgressBar::operator++(){
    *this += 1;
    return *this;
}

ConsoleProgressBar ConsoleProgressBar::operator++(int){
    ConsoleProgressBar tmp(*this);
    operator++();
    return tmp;
}

void ConsoleProgressBar::complete(){
    step(mExpectedCount - mCount);
}

void ConsoleProgressBar::start(){
    draw();
}


