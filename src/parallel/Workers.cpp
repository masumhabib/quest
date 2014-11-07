/* 
 * File:   Workers.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 16, 2014, 6:04 PM
 */

#include "parallel/Workers.h"


namespace qmicad{namespace parallel{

Workers::Workers():mMasterId(0), mWorkers(mworld)
{
    init();
}

Workers::Workers(const vector<string> &argv):mMasterId(0), mWorkers(mworld)
{
}

Workers::Workers(int argc, char** argv):mMasterId(0), menv(argc, argv), mWorkers(mworld)
{
    init();
}

Workers::Workers( const communicator &workers):mWorkers(workers), mMasterId(0)
{
    init();
}

void Workers::init(){
    mMyCpuId = mWorkers.rank();
    mNcpu = mWorkers.size();
    mIAmMaster = (mMasterId == mMyCpuId);

    stds::vout.printersId(mMasterId);        
    stds::vout.myId(mMyCpuId);    
}

void Workers::assignCpus(long& myStart, long& myEnd, long& myN, long N) const{
    
    long quo = N/mNcpu;
    long rem = N%mNcpu;
    myStart = mMyCpuId*quo + (mMyCpuId < rem ? mMyCpuId:rem);
    myEnd = myStart + quo + (mMyCpuId < rem ? 1:0) - 1;        
    myN = myEnd - myStart + 1;
}
}}



