/* 
 * File:   Workers.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 16, 2014, 6:04 PM
 */

#include "parallel/Workers.h"
#include "python/boostpython.hpp"


namespace qmicad{
namespace python{

Workers::Workers( const communicator &workers):mWorkers(workers), mMasterId(0)
{
    mMyCpuId = mWorkers.rank();
    mNcpu = mWorkers.size();
    mIAmMaster = (mMasterId == mMyCpuId);

    stds::vout.printersId(mMasterId);        
    stds::vout.myId(mMyCpuId);        
}

void Workers::assignCpus(long& myStart, long& myEnd, long& myN, long N){
    
    long quo = N/mNcpu;
    long rem = N%mNcpu;
    myStart = mMyCpuId*quo + (mMyCpuId < rem ? mMyCpuId:rem);
    myEnd = myStart + quo + (mMyCpuId < rem ? 1:0) - 1;        
    myN = myEnd - myStart + 1;
}

/**
 * MPI communicator wrapper.
 */
void export_Workers(){
    using namespace utils;
    
   class_<Workers, shared_ptr<Workers>, noncopyable>("Workers", 
          init<const communicator&>())
        .def("MyId", &Workers::MyId)
        .def("MasterId", &Workers::MasterId)
        .def("N", &Workers::N)
        .def("AmIMaster", &Workers::AmIMaster)
        .def("IAmMaster", &Workers::IAmMaster)
    ;     
}

}
}

