/* 
 * File:   Workers.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 16, 2014, 6:03 PM
 */

#ifndef WORKERS_H
#define	WORKERS_H

#include "utils/vout.h"
#include <boost/mpi.hpp>

namespace utils{
using namespace boost::mpi;

class Workers {
public:
    Workers(const communicator &workers):
            mWorkers(workers), mMasterId(0)
    {
        mMyCpuId = mWorkers.rank();
        mNcpu = mWorkers.size();
        mIAmMaster = (mMasterId == mMyCpuId);

        stds::vout.printersId(mMasterId);        
        stds::vout.myId(mMyCpuId);        
    }    

    int     MyId()          const {return mMyCpuId; };
    int     MasterId()      const { return mMasterId; };
    int     N()             const { return mNcpu; };
    bool    AmIMaster()     const { return mIAmMaster; };
    bool    IAmMaster()     const { return mIAmMaster; };
    
    const communicator& Comm() const {return mWorkers; };
    
protected:
    const communicator      &mWorkers;      //!< MPI worker processes 
    int                     mMyCpuId;       //!< This process ID
    int                     mNcpu;          //!< Total number of processes
    bool                    mIAmMaster;     //!< Master process
    const int               mMasterId;      //!< Master ID

};


}
#endif	/* WORKERS_H */

