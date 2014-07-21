/* 
 * File:   NegfEloop.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "negf/CohRgfLoop.h"
#include "negf/NegfResult.h"
#include "python/boostpython.hpp"

namespace qmicad{
namespace negf{

CohRgfLoop::CohRgfLoop(const Workers &workers, uint nb, double kT, dcmplx ieta, 
        bool orthogonal, string newprefix): Printable(newprefix), 
        mrgf(nb, kT, ieta, orthogonal, " " + newprefix), mbar("  NEGF: "),
        mWorkers(workers) , mnb(nb)       
{
    
}

void CohRgfLoop::E(const vec &E){
    mE = E;
    mbar.expectedCount(npoints());
}
 
void CohRgfLoop::k(const mat &k){
    mk = k;
    mbar.expectedCount(npoints());
}

void CohRgfLoop::mu(double muD, double muS){
    mrgf.mu(muD, muS);
}

void CohRgfLoop::run(){
    
    long n = npoints();
    long nE = mE.is_empty()?1:mE.n_rows;
    long nk = mk.is_empty()?1:mk.n_rows;
    
    
    long myStart, myEnd, myN;
    mWorkers.assignCpus(myStart, myEnd, myN, n);
    for(long it = myStart; it < myEnd; ++it){
        long ik, iE;
        field<shared_ptr<cxmat> > H0(mnb);
        field<shared_ptr<cxmat> > S0(mnb);
        field<shared_ptr<cxmat> > Hl(mnb+1);
        field<shared_ptr<cxmat> > Sl(mnb+1);
        if (!mk.is_empty()){
            ik = it/nE;
            iE = it%nE;
            
            vec k = mk.row(ik);
            // generate H0(k) = sum H0(n)*exp(i k.rn)
            for(int ib = 0; ib < mnb; ++ib){
                shared_ptr<cxmat> H0k = make_shared(new cxmat(mH0(0,ib)->n_rows, mH0(0, ib)->n_cols, fill::zeros));
                shared_ptr<cxmat> S0k = make_shared(new cxmat(mS0(0,ib)->n_rows, mS0(0, ib)->n_cols, fill::zeros));
                for(int in = 0; in < mH0.n_rows; ++in){
                    cxmat &Hk = *H0k;
                    double th = dot(k, mpv0(in));
                    Hk = Hk + (*mH0(in, ib))*exp(i*th);
                    cxmat &Sk = *S0k;
                    Sk = Sk + (*mS0(in, ib))*exp(i*th);
                }
                H0(ib) = H0k;
                S0(ib) = S0k;
            }
            
            
            
        }else{
            ik = 0;
            iE = it;
        }
        
        // set energy
        mrgf.E(mE[iE]);
        
        // run simulation
        prepare();
        preCompute();
        compute();
        postCompute();
        collect();
    }
    
}

void CohRgfLoop::prepare() {
    mWorkers.Comm().barrier();
    mbar.start();
}

void CohRgfLoop::preCompute(int il){
    double E = mE(il);
    //mnegf = shared_ptr<CohRgfa>(new CohRgfa(mnp, E));
};

void CohRgfLoop::compute(int il){
    negf_result r;                       // result as a function of energy
    r.E = mE(il);                  
    // Transmission
    if(mTE.isEnabled()){
        r.M = mrgf.TEop(mTE.N);  // M => T(E)
        mThisTE.push_back(r);  
    }
    // Current
    for (int it = 0; it < mIop.size(); ++it){
        r.M = mrgf.Iop(mIop[it].N,  mIop[it].ib, mIop[it].jb); 
        mThisIop[it].push_back(r);           // ThisIop[it] => vector of Iop()
    }
    // Density of States
    if(mDOS.isEnabled()){
        r.M = mrgf.DOSop(mDOS.N);  // M => DOS(E)
        mThisDOS.push_back(r);  
    }
    // Electron density
    if(mn.isEnabled()){
        r.M = mrgf.nop(mn.N);  // M => n(E)
        mThisn.push_back(r);  
    }

}

void CohRgfLoop::postCompute(int il){
    mrgf.reset();          // free up memory
    ++mbar;                 // Show feedback
}

void CohRgfLoop::collect(){
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

void CohRgfLoop::gather(vec_result &thisR, NegfResultList &all){

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

long CohRgfLoop::npoints(){
    long n = (mE.is_empty()?1:mE.n_rows) * (mk.is_empty()?1:mk.n_rows);
    return n;
}

void CohRgfLoop::save(string fileName){
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

void CohRgfLoop::enableTE(uint N){
    mTE.tag = "TRANSMISSION";
    mTE.N = N;
}

void CohRgfLoop::enableI(uint N, uint ib, uint jb){
    stringstream out;
    out << "CURRENT";
    mIop.push_back(NegfResultList(out.str(), N, ib, jb));
    mThisIop.push_back(vec_result());
}

void CohRgfLoop::enableDOS(uint N){
    mDOS.tag = "DOS";
    mDOS.N = N;
}

void CohRgfLoop::enablen(uint N){
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
    class_<CohRgfLoop, shared_ptr<CohRgfLoop> >("NegfEloop", 
            init<VecGrid&, /*const CohRgfaParams&,*/ const Workers&, 
            optional<bool> >())
        .def("run", &CohRgfLoop::run)
        .def("save", &CohRgfLoop::save)
        .def("enableTE", &CohRgfLoop::enableTE, NegfEloop_enableTE())
        .def("enableI", &CohRgfLoop::enableI, NegfEloop_enableI())
        .def("enableDOS", &CohRgfLoop::enableDOS, NegfEloop_enableDOS())
        .def("enablen", &CohRgfLoop::enablen, NegfEloop_enablen())
    ;
}

}
}



