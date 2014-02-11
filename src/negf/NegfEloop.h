/* 
 * File:   NegfEloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#ifndef NEGFELOOP_H
#define	NEGFELOOP_H

#include "../utils/std.hpp"
#include "../utils/serialize.hpp"
#include "../maths/fermi.hpp"
#include "../parallel/parloop.h"


#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/smart_ptr.hpp>

#include "CohRgfa.h"

namespace qmicad{
using namespace utils::stds;
using utils::ParLoop;
using utils::VecGrid;
using boost::shared_ptr;
namespace mpi = boost::mpi;

/*
 * NegfCalculations specifies what to calculate.
 */
class NegfEloop: public ParLoop<double>{
public:
    NegfEloop(const VecGrid &E, const NegfParams &np, const mpi::communicator &workers):
            ParLoop<double>(E, workers),TE(1), I1(2), mnp(np)
    {
        mTEn = 1;
        mCalcType = TE;
        mprog   = 0;
        mprogmx = 80;
    }
    
    void            enableTE(uint N = 1){mCalcType |= TE; mTEn = N; };
    virtual void    saveTE(string FileName);
    
    
protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    
    virtual void    computeTE(uint N = 1);
    virtual void    collectTE();
    
    virtual void    stepCompleted();
    virtual void    gather(vector<cxmat> &Mproc, list<cxmat> &M);

public:
    const uint           TE;
    const uint           I1;
    
protected:
    const NegfParams      &mnp;     // Negf parameters
    shared_ptr<CohRgfa>      mnegf;    // current Negf calculator
    double                mEi;      // current energy
    
    // output
    uint                 mCalcType; // Which calculations we need to perform
    // TE
    vector<cxmat>        mTEproc;  // transmission list for local process
    list<cxmat>          mTE;      // transmission list for all processes
    uint                 mTEn;     // size of TE output matrix
    
    
        
    // user feedback
    int                   mprog;    // progress
    int                   mprogmx;  // maximum progress
    
};
}
#endif	/* NEGFELOOP_H */

