/* 
 * File:   NegfEloop.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "NegfEloop.h"
#include "NegfResult.h"
#include "../python/boostpython.hpp"

namespace qmicad{
namespace negf{

NegfEloop::NegfEloop(const VecGrid &E, const CohRgfaParams &np, 
        const Workers &workers, bool isAscii): ParLoop(workers, E.N()), 
        mnp(np), mE(E), mIsAscii(isAscii), mbar("  NEGF: ",  E.N()/workers.N()+1)
{   
}
    

void NegfEloop::prepare() {
    mWorkers.Comm().barrier();
    mbar.start();
}

void NegfEloop::preCompute(int il){
    double E = mE(il);
    mnegf = shared_ptr<CohRgfa>(new CohRgfa(mnp, E));
};

void NegfEloop::compute(int il){
    negf_result r;                       // result as a function of energy
    r.E = mE(il);                  
    // Transmission
    if(mTE.isEnabled()){
        r.M = mnegf->TEop(mTE.N);  // M => T(E)
        mThisTE.push_back(r);  
    }
    // Current
    for (int it = 0; it < mIop.size(); ++it){
        r.M = mnegf->Iop(mIop[it].N,  mIop[it].ib, mIop[it].jb); 
        mThisIop[it].push_back(r);           // ThisIop[it] => vector of Iop()
    }
    // Density of States
    if(mDOS.isEnabled()){
        r.M = mnegf->DOSop(mDOS.N);  // M => DOS(E)
        mThisDOS.push_back(r);  
    }
    // Electron density
    if(mn.isEnabled()){
        r.M = mnegf->nop(mn.N);  // M => n(E)
        mThisn.push_back(r);  
    }

}

void NegfEloop::postCompute(int il){
    mnegf.reset();          // free up memory
    ++mbar;                 // Show feedback
}

void NegfEloop::collect(){
    // Update the progress bar.
    mWorkers.Comm().barrier();
    mbar.complete();
    
    // Gather transmission
    if(mTE.isEnabled()){
        gather(mThisTE, mTE);
    }

    // Gather current
    for (int it = 0; it < mIop.size(); ++it){
        gather(mThisIop[it], mIop[it]);
    }

    // Gather Density of States
    if(mDOS.isEnabled()){
        gather(mThisDOS, mDOS);
    }

    // Gather electron density
    if(mn.isEnabled()){
        gather(mThisn, mn);
    }

}

void NegfEloop::gather(vec_result &thisR, NegfResultList &all){

    if(!mWorkers.IAmMaster()){    
        // slaves send their local data
        mpi::gather(mWorkers.Comm(), thisR, mWorkers.MasterId());

    // The master collects data        
    }else{
        // Collect T(E)
        vector<vector<negf_result> >gatheredR(mN);
        mpi::gather(mWorkers.Comm(), thisR, gatheredR, mWorkers.MasterId());
        
        // merge and store results on mTE list.
        vector<vector<negf_result> >::iterator it;
        for (it = gatheredR.begin(); it != gatheredR.end(); ++it){
            all.R.insert(all.R.end(), it->begin(), it->end());
        }
        // Sort the results based on energy.
        all.sort();
    }       
}

void NegfEloop::save(string fileName){
    if(mWorkers.IAmMaster()){
        // save to a file
        ofstream out;
        if(mIsAscii){
            out.open(fileName.c_str(), ostream::binary|ios::app);
        }else{
            out.open(fileName.c_str(), ios::app);
        }
        if (!out.is_open()){
            throw ios_base::failure(" NegfResult::saveTE(): Failed to open file " 
                    + fileName + ".");
        }
        
        // Transmission 
        if (mTE.isEnabled()){
            mTE.save(out);
        }
        
        // Current
        for (int it = 0; it < mIop.size(); ++it){
            mIop[it].save(out);
        }

        // Density of States
        if (mDOS.isEnabled()){
            mDOS.save(out);
        }

        // Electron density
        if (mn.isEnabled()){
            mn.save(out);
        }

    }
}

void NegfEloop::enableTE(uint N){
    mTE.tag = "TRANSMISSION";
    mTE.N = N;
}

void NegfEloop::enableI(uint N, uint ib, uint jb){
    stringstream out;
    out << "CURRENT";
    mIop.push_back(NegfResultList(out.str(), N, ib, jb));
    mThisIop.push_back(vec_result());
}

void NegfEloop::enableDOS(uint N){
    mDOS.tag = "DOS";
    mDOS.N = N;
}

void NegfEloop::enablen(uint N){
    mn.tag = "n";
    mn.N = N;
}

}
}

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace negf;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableTE, enableTE, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableI, enableI, 0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enableDOS, enableDOS, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(NegfEloop_enablen, enablen, 0, 1)
void export_NegfEloop(){
    class_<NegfEloop, shared_ptr<NegfEloop> >("NegfEloop", 
            init<VecGrid&, const CohRgfaParams&, const Workers&, 
            optional<bool> >())
        .def("run", &NegfEloop::run)
        .def("save", &NegfEloop::save)
        .def("enableTE", &NegfEloop::enableTE, NegfEloop_enableTE())
        .def("enableI", &NegfEloop::enableI, NegfEloop_enableI())
        .def("enableDOS", &NegfEloop::enableDOS, NegfEloop_enableDOS())
        .def("enablen", &NegfEloop::enablen, NegfEloop_enablen())
    ;
}

}
}



