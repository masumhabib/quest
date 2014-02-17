/* 
 * File:   parloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 2:41 PM
 */

#ifndef PARFOR_H
#define	PARFOR_H

#include "../grid/grid.hpp"
#include "../utils/vout.h"
#include "Workers.h"


namespace utils{
/*
 * Parallel loop
 */
class ParLoop{
public: 
    ParLoop(const Workers &workers, int N);
    virtual void assign(int N);
    virtual void run();
    
protected:
    virtual void prepare(){};
    virtual void preCompute(int il){};
    virtual void compute(int il){};    
    virtual void postCompute(int il){};
    virtual void collect(){};
    
    
protected:
    int                 mN;             //!< Number of grid points.
    int                 mMyN;           //!< Number of grid points to be calculated by this process.
    const Workers       &mWorkers;      //!< MPI worker processes.
    
    int                 mMyStart;       //!< Start point for this CPU.
    int                 mMyEnd;         //!< End point for this CPU.
};

}

#endif	/* PARFOR_H */

