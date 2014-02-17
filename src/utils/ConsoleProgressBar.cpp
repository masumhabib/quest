/* 
 * File:   ConsoleProgressBar.cpp
 * Author: masum
 * 
 * Created on April 8, 2013, 11:35 AM
 */

#include "ConsoleProgressBar.h"    

namespace utils{
    

ConsoleProgressBar::ConsoleProgressBar(string prefix, unsigned long expectedCount, 
    unsigned long count, bool livePercent):
    mExpectedCount(expectedCount),
    mCount(count < expectedCount ? count:expectedCount),
    mPrefix(prefix), mLivePercent(livePercent)
{
    mSuffix = "%";
    mSpace = ' ';
    mDot   = '=';
    mEnds  = '|';

    mnDots = convertCountToDots();
}

ConsoleProgressBar::~ConsoleProgressBar() {
}

void ConsoleProgressBar::draw(){
    // Updates live progress percentage.
    if(mLivePercent){
        string bar;
        stringstream ssbar;

        bar.insert(bar.length(), mnDots, mDot);
        bar.insert(bar.length(), mnTotalDots - mnDots, mSpace);

        ssbar << mPrefix << "[";
        ssbar.width(3);
        ssbar << int(mnDots*mPercentPerDot);
        ssbar << mSuffix;
        ssbar << "] " << mEnds;
        ssbar << bar << mEnds;

        vout << ssbar.str() << endl;
        vout << "\033[F";
    // Does not show the live progress percentage.
    }else{
        if (mnDots == 0){
            vout << vnormal << mPrefix << mEnds;
        }else if(mnDots == mnTotalDots){
            vout << vnormal << mEnds;
        }else{
            vout << vnormal << mDot;
        }
        
        
    }
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

}

