/* 
 * File:   NegfEloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#ifndef NEGFELOOP_H
#define	NEGFELOOP_H

#include "../utils/serialize.hpp"
#include "../utils/fermi.hpp"
#include "../parallel/parloop.h"

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/smart_ptr.hpp>

#include "NEGF.h"

namespace qmicad{
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
    
    void            enableTE(uword N = 1){mCalcType |= TE; mTEn = N; };
    virtual void    saveTE(string FileName);
    
    
protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    
    virtual void    computeTE(uword N = 1);
    virtual void    collectTE();
    
    virtual void    stepCompleted();
    virtual void    gather(vector<cx_mat> &Mproc, list<cx_mat> &M);

public:
    const uword           TE;
    const uword           I1;
    
protected:
    const NegfParams      &mnp;     // Negf parameters
    shared_ptr<Negf>      mnegf;    // current Negf calculator
    double                mEi;      // current energy
    
    // output
    uword                 mCalcType; // Which calculations we need to perform
    // TE
    vector<cx_mat>        mTEproc;  // transmission list for local process
    list<cx_mat>          mTE;      // transmission list for all processes
    uword                 mTEn;     // size of TE output matrix
    
    
        
    // user feedback
    int                   mprog;    // progress
    int                   mprogmx;  // maximum progress
    
};
}
#endif	/* NEGFELOOP_H */

