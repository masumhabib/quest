/* 
 * File:   NegfEloop.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#ifndef NEGFELOOP_H
#define	NEGFELOOP_H

#include "negf/CohRgfa.h"
#include "negf/NegfResult.h"

#include "utils/ConsoleProgressBar.h"
#include "utils/std.hpp"
#include "utils/vout.h"
#include "utils/serialize.hpp"
#include "maths/fermi.hpp"
#include "parallel/parloop.h"


#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/smart_ptr.hpp>

namespace qmicad{
namespace negf{

using namespace utils::stds;
using utils::ParLoop;
using utils::VecGrid;
using boost::shared_ptr;
namespace mpi = boost::mpi;

/*
 * NegfCalculations specifies what to calculate.
 */
class NegfEloop: public ParLoop{
    typedef vector<negf_result> vec_result;
public:
    NegfEloop(const VecGrid &E, const CohRgfaParams &np, const Workers &workers, 
              bool isText = true);
    
    void            enableTE(uint N = 1);   
    void            enableI(uint N = 1, uint ib = 0, uint jb = 0);
    void            enableDOS(uint N = 1);
    void            enablen(uint N = 1, int ib = -1); //!< Electron density.
    void            enablep(uint N = 1, int ib = -1); //!< Hole density.
    
    virtual void    save(string fileName);
    
protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    virtual void    gather(vec_result &thisR, NegfResultList &all);

public:
    
protected:
    const CohRgfaParams   &mnp;         //!< Negf parameters.
    shared_ptr<CohRgfa>   mnegf;        //!< Current Negf calculator.
    VecGrid               mE;           //!< Energy grid.
    
    // TEop: Transmission operator
    vec_result            mThisTE;      //!< Transmission list for local process
    NegfResultList        mTE;          //!< Transmission list for all processes

    // Iop, Current operator for block # i.
    vector<vec_result>  mThisIop;       //!< For local process.
    vector<NegfResultList> mIop;        //!< Collection of all processes.
    
    // DOSop: Density of States operator
    vec_result            mThisDOS;     //!< DOS list for local process
    NegfResultList        mDOS;         //!< DOS list for all processes
    
    // Electron density operator
    vector<vec_result>    mThisnOp;     //!< Density list for local process
    vector<NegfResultList>mnOp;         //!< Density list for all processes

    // Hole density operator
    vector<vec_result>    mThispOp;     //!< Density list for local process
    vector<NegfResultList>mpOp;         //!< Density list for all processes    
    
    bool                  mIsText;     //!< Save as ascii/binary?

    // user feedback
    ConsoleProgressBar    mbar;         //!< Shows a nice progress bar.
    
};
}
}
#endif	/* NEGFELOOP_H */

