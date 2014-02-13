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
#include "NegfResult.h"


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
    NegfEloop(const VecGrid &E, const NegfParams &np, 
              const mpi::communicator &workers, bool saveAscii = true);
    
    void            enableTE(uint N = 1);
    void            disableTE();    
    virtual void    save(string fileName);
    
protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    virtual void    stepCompleted();
    virtual void    gather(vector<negfresult> &thisR, NegfResultList &all);

public:
    
protected:
    const NegfParams      &mnp;         //!< Negf parameters.
    shared_ptr<CohRgfa>   mnegf;        //!< Current Negf calculator.
    
    // TE
    vector<negfresult>    mThisTE;      //!< Transmission list for local process
    NegfResultList        mTE;          //!< Transmission list for all processes
    
    // user feedback
    int                  mprog;         //!< Calculation progress.
    int                  mprogmx;       //!< Maximum progress.
    
};
}
#endif	/* NEGFELOOP_H */

