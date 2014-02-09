/* 
 * File:   parloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 2:41 PM
 */

#ifndef PARFOR_H
#define	PARFOR_H

#include <boost/mpi.hpp>
#include "../grid/grid.hpp"

namespace utils{
namespace mpi = boost::mpi;

/*
 * Parallel loop
 */
template<class T>
class ParLoop{
public: 
    ParLoop(const Grid1D<T> &L, const mpi::communicator &workers):
            mL(L), mWorkers(workers), mMasterId(0)
    {
        mMyCpuId = mWorkers.rank();
        mNcpu = mWorkers.size();
        mIAmMaster = (mMasterId == mMyCpuId);
        mN = mL.N();
        
        int quo = mN/mNcpu;
        int rem = mN%mNcpu;
        mMyStart = mMyCpuId*quo + (mMyCpuId < rem ? mMyCpuId:rem);
        mMyEnd = mMyStart + quo + (mMyCpuId < rem ? 1:0) - 1;
        
    }    

    virtual void run(){
        prepare();
        for(int il = mMyStart; il <= mMyEnd; ++il){
            preCompute(il);
            compute(il);
            postCompute(il);
        }
        collect();
    }
protected:
    virtual void prepare(){};
    virtual void preCompute(int il){};
    virtual void compute(int il){};    
    virtual void postCompute(int il){};
    virtual void collect(){};
    
    
protected:
    Grid1D<T>           mL;             // loop variable
    int                 mN;             // number of grid points
    const mpi::communicator &mWorkers;       // MPI worker processes 
    int                 mMyCpuId;       // This process ID
    int                 mNcpu;          // Total number of processes
    bool                mIAmMaster;     // Master process
    const int           mMasterId;      // Master ID
    
    int                 mMyStart;       // Start point for this CPU
    int                 mMyEnd;         // End point for this CPU
};

}

#endif	/* PARFOR_H */

