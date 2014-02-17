/* 
 * File:   NegfEloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#ifndef NEGFELOOP_H
#define	NEGFELOOP_H

#include "NegfResult.h"

#include "../utils/ConsoleProgressBar.h"
#include "../utils/std.hpp"
#include "../utils/vout.h"
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
class NegfEloop: public ParLoop{
public:
    NegfEloop(const VecGrid &E, const NegfParams &np, const Workers &workers, 
              bool saveAscii = true);
    
    void            enableTE(uint N = 1);
    void            disableTE();    
    void            enableI1(uint N = 1);
    void            enableI1sx(uint N = 1);
    void            enableI1sy(uint N = 1);
    void            enableI1sz(uint N = 1);
    
    virtual void    save(string fileName);
    
protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    virtual void    gather(vector<negfresult> &thisR, NegfResultList &all);

public:
    
protected:
    const NegfParams      &mnp;         //!< Negf parameters.
    shared_ptr<CohRgfa>   mnegf;        //!< Current Negf calculator.
    VecGrid               mE;           //!< Energy grid.
    
    // TE
    vector<negfresult>    mThisTE;      //!< Transmission list for local process
    NegfResultList        mTE;          //!< Transmission list for all processes

    // I1op                             //!< Current operator for terminal # 1.
    vector<negfresult>    mThisI1op;    //!< For local process.
    NegfResultList        mI1op;        //!< Collection of all processes.
    uint                  mI1N;         //!< Calculate charge current?
    uint                  mI1sxN;       //!< Calculate spin x current?                   
    uint                  mI1syN;       //!< Calculate spin y current?
    uint                  mI1szN;       //!< Calculate spin z current?
    
    // user feedback
//    int                   mprog;         //!< Calculation progress.
//    int                   mprogmx;       //!< Maximum progress.
    ConsoleProgressBar    mbar;          //!< Shows a nice progress bar.
    
};
}
#endif	/* NEGFELOOP_H */

