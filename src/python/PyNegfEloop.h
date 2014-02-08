/* 
 * File:   PyNegfEloop.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 8, 2014, 12:15 AM
 * 
 * Python wrapper for NegfEloop class
 * 
 */

#ifndef PYNEGFELOOP_H
#define	PYNEGFELOOP_H

#include "../negf/NegfEloop.h"

namespace qmicad{namespace python{
/**
 * Python wrapper for NegfEloop class.
 */
class PyNegfEloop:public NegfEloop{
public:
   
private:

public:
    PyNegfEloop(const VecGrid &E, const NegfParams &np, const mpi::communicator &workers):
        NegfEloop(E, np, workers){
    };
    
    virtual void    prepare(){ NegfEloop::prepare(); };
    virtual void    preCompute(int il) { NegfEloop::preCompute(il); };
    virtual void    compute(int il) { NegfEloop::compute(il); };  
    virtual void    postCompute(int il) { NegfEloop::postCompute(il); };
    virtual void    collect() { NegfEloop::collect(); };
    
    virtual void    computeTE(uword N = 1) {NegfEloop::computeTE(N); };
    virtual void    collectTE() {NegfEloop::collectTE(); };
    
    virtual void    stepCompleted() {NegfEloop::stepCompleted(); };
    
};
}}
#endif	/* PYNEGFELOOP_H */

