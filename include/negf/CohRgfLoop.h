/* 
 * File:   CohRgfLoop.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#ifndef COHRGFLOOP_H
#define	COHRGFLOOP_H

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
using boost::make_shared;
namespace mpi = boost::mpi;

/*
 * NegfCalculations specifies what to calculate.
 */
class CohRgfLoop: public Printable{
    typedef vector<negf_result> vec_result;
public:
    CohRgfLoop(const Workers &workers, uint nb = 5, double kT = 0.0259, 
        dcmplx ieta = dcmplx(0,1E-3), bool orthogonal = true, string newprefix = ""); 
    //const VecGrid &E, /*const CohRgfaParams &np,*/ const Workers &workers,  bool isAscii = true);

    void            E(const vec &E);
    void            k(const mat &k);
    void            mu(double muD = 0.0, double muS = 0.0);
    
    void            enableTE(uint N = 1);   
    void            enableI(uint N = 1, uint ib = 0, uint jb = 0);
    void            enableDOS(uint N = 1);
    void            enablen(uint N = 1);
    
    void            run();
    virtual void    save(string fileName);
    
private:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    virtual void    gather(vec_result &thisR, NegfResultList &all);
    
    long            npoints();

public:
    
public:
    const Workers         &mWorkers;    //!< MPI worker processes.
    CohRgfa               mrgf;         //!< Current Negf calculator.
    vec                   mE;           //!< Energy grid.
    mat                   mk;          //!< Wave vector.
    int                   mnb;
    // Hamiltonian , overlap and potential
    field<shared_ptr<cxmat> >mH0;// Diagonal blocks of Hamiltonian: H0(i) = [H]_i,i
                                 // H0(0) is on the left contact and H0(N+1) is 
                                 // on the right contact.
    field<vec>              mpv0;// Position vectors of block H0(i)
    field<shared_ptr<cxmat> >mS0;// Diagonal blocks of overlap matrix: S0(i) = [S]_i,i
                                 // H0(0) is on the left contact and H0(N+1) is 
                                 // on the right contact.    
    field<shared_ptr<cxmat> >mHl;// Lower diagonal blocks of Hamiltonian: 
                                 // Hl(i) = [H]_i,i-1. Hl(0) = [H]_0,-1 is the 
                                 // hopping between two left contact blocks.
                                 // H(1) = [H]_1,0 is the hopping 
                                 // between block # 1 and left contact. Hl(N+1) is
                                 // hopping between right contact and block # N.
                                 // Hl(N+2) is hopping between two right contact
                                 // blocks.
    field<vec>              mpvl;// Position vectors of block Hl(i)
    field<shared_ptr<cxmat> >mSl;// Lower diagonal blocks of overlap matrix: 
                                 // Sl(i) = [S]_i,i-1. Sl(0) = [S]_0,-1 is the 
                                 // overlap between two left contact blocks.
                                 // S(1) = [S]_1,0 is the overlap between 
                                 // block # 1 and left contact. Hl(N+1) is
                                 // overlap between right contact and block # N.
                                 // Sl(N+2) is hopping between two right contact
                                 // blocks.
    field<shared_ptr<vec> >   mV;// Electrostatic potential of all the orbitals
                                 // for the entire device: from block#0 
                                 // to block#N+1.
    
    
    // TEop: Transmission operator
    vec_result            mThisTE;      //!< Transmission list for local process.
    NegfResultList        mTE;          //!< Transmission list for all processes.

    // Iop, Current operator for block # i.
    vector<vec_result>  mThisIop;       //!< For local process.
    vector<NegfResultList> mIop;        //!< Collection of all processes.
    
    // DOSop: Density of States operator
    vec_result            mThisDOS;     //!< DOS list for local process.
    NegfResultList        mDOS;         //!< DOS list for all processes.
    
    // Electron density operator
    vec_result            mThisn;     //!< Density list for local process.
    NegfResultList        mn;         //!< Density list for all processes.
    
    bool                  mIsAscii;     //!< Save as ascii/binary?

    // user feedback
    ConsoleProgressBar    mbar;         //!< Shows a nice progress bar.
    
};
}
}
#endif	/* COHRGFLOOP_H */

