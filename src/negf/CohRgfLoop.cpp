/* 
 * File:   NegfEloop.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "negf/CohRgfLoop.h"
#include "negf/RgfResult.h"
#include "python/boostpython.hpp"

namespace qmicad{
namespace negf{

CohRgfLoop::CohRgfLoop(const Workers &workers, uint nb, double kT, dcmplx ieta, 
        bool orthogonal, string newprefix): Printable(newprefix), 
        mrgf(nb, kT, ieta, orthogonal, " " + newprefix), mbar("  NEGF: "),
        mWorkers(workers) , mnb(nb), morthogonal(morthogonal)
{
    
}

void CohRgfLoop::E(const vec &E){
    if (E.is_empty()){
        throw runtime_error("In CohRgfLoop::E(), E cannot be empty.");
    }
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
    
    prepare();

    long n = npoints();
    long nE = mE.n_rows;
    long nk = mk.n_rows;
    
    long myStart, myEnd, myN;
    // Assign E and k points to CPUs 
    mWorkers.assignCpus(myStart, myEnd, myN, n);
    // Loop over problem assigned to this CPU.
    for(long it = myStart; it < myEnd; ++it){
        
        long ik, iE;
        // Hamiltonian and overlap matrices. 
        field<shared_ptr<cxmat> > H0(mnb);
        field<shared_ptr<cxmat> > S0(mnb);
        field<shared_ptr<cxmat> > Hl(mnb+1);
        field<shared_ptr<cxmat> > Sl(mnb+1);

        if (nk != 0){ // Do a k-loop
            //@TODO: Needs memory optimization.
            ik = it/nE;
            iE = it%nE;
            
            vec k = mk.row(ik);
            
            // Calculate block diagonal matrices.
            for(int ib = 0; ib < mnb; ++ib){
                shared_ptr<cxmat> H0k = make_shared<cxmat>(mH0(0,ib)->n_rows, mH0(0, ib)->n_cols, fill::zeros);
                shared_ptr<cxmat> S0k = make_shared<cxmat>(mS0(0,ib)->n_rows, mS0(0, ib)->n_cols, fill::zeros);
                
                double th;
                dcmplx expith;
                for(int in = 0; in < mH0.n_rows; ++in){
                    th = dot(k, mpv0(in));
                    expith = exp(i*th);
                    (*H0k) = (*H0k) + (*mH0(in, ib))*expith;                    
                    (*S0k) = (*S0k) + (*mS0(in, ib))*expith;
                }
                H0(ib) = H0k;
                S0(ib) = S0k;
            }
            // Calculate lower block diagonals
            for(int ib = 0; ib <= mnb; ++ib){
                shared_ptr<cxmat> Hlk = make_shared<cxmat>(mHl(0,ib)->n_rows, mHl(0, ib)->n_cols, fill::zeros);
                shared_ptr<cxmat> Slk = make_shared<cxmat>(mSl(0,ib)->n_rows, mSl(0, ib)->n_cols, fill::zeros);
                
                double th;
                dcmplx expith;
                for(int in = 0; in < mHl.n_rows; ++in){
                    th = dot(k, mpvl(in));
                    expith = exp(i*th);
                    (*Hlk) = (*Hlk) + (*mHl(in, ib))*expith;
                    if (!morthogonal){
                        (*Slk) = (*Slk) + (*mSl(in, ib))*expith;
                    }
                }
                Hl(ib) = Hlk;
                if (!morthogonal){
                    Sl(ib) = Slk;
                }
            }
        }else{ // Do only E loop
            ik = 0;
            iE = it;
            H0 = mH0;
            S0 = mS0;
            Hl = mHl;
            Sl = mSl;   
        }
        
        // set E, H and S.
        mrgf.E(mE[iE]);
        mrgf.H(H0, Hl);
        mrgf.S(S0, Sl);
        mrgf.V(mV);
        
        // run simulation step.
        compute();
        ++mbar;                // Show feedback        
    }
    
    collect();
    
}

void CohRgfLoop::prepare() {
    mWorkers.Comm().barrier();
    mbar.start();
}

void CohRgfLoop::compute(){
    cxmat r;                       // result as a function of energy
    
    // Transmission
    if(mTE.isEnabled()){
        r = mrgf.TEop(mTE.N);  // M => T(E)
        mThisTE.push_back(r);  
    }
    // Current
    for (int it = 0; it < mIop.size(); ++it){
        r = mrgf.Iop(mIop[it].N,  mIop[it].ib, mIop[it].jb); 
        mThisIop[it].push_back(r);           // ThisIop[it] => vector of Iop()
    }
    // Density of States
    if(mDOS.isEnabled()){
        r = mrgf.DOSop(mDOS.N);  // M => DOS(E)
        mThisDOS.push_back(r);  
    }
    // Non-equilibrium electron density
    for (int it = 0; it < mnOp.size(); ++it){
        r = mrgf.nOp(mnOp[it].N,  mnOp[it].ib); 
        mThisnOp[it].push_back(r);           // Thisnop[it] => vector of nop()
    }
    // Non-equilibrium hole density
    for (int it = 0; it < mpOp.size(); ++it){
        r = mrgf.pOp(mpOp[it].N,  mpOp[it].ib); 
        mThispOp[it].push_back(r);        
    }    

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

    // Gather equilibrium electron density
    for (int it = 0; it < mnOp.size(); ++it){
        gather(mThisnOp[it], mnOp[it]);
    }    

    // Gather Non-equilibrium electron density
    for (int it = 0; it < mpOp.size(); ++it){
        gather(mThisnOp[it], mpOp[it]);
    }    

}

void CohRgfLoop::gather(cxmat_vec &thisR, RgfResult &all){

    if(!mWorkers.IAmMaster()){    
        // slaves send their local data
        mpi::gather(mWorkers.Comm(), thisR, mWorkers.MasterId());

    // The master collects data        
    }else{
        // Collect T(E)
        vector<cxmat_vec>gatheredR(mWorkers.N());
        mpi::gather(mWorkers.Comm(), thisR, gatheredR, mWorkers.MasterId());
        
        // merge and store results on mTE list.
        vector<cxmat_vec>::iterator it;
        for (it = gatheredR.begin(); it != gatheredR.end(); ++it){
            all.R.insert(all.R.end(), it->begin(), it->end());
        }
    }       
}

long CohRgfLoop::npoints(){
    long n = (mE.is_empty()?1:mE.n_rows) * (mk.is_empty()?1:mk.n_rows);
    return n;
}

void CohRgfLoop::save(string fileName, bool isText){
    if(mWorkers.IAmMaster()){
        // save to a file
        
        if(isText){ // ASCII format
            ofstream out;
            out.open(fileName.c_str(), ios::app);
            if (!out.is_open()){
                throw ios_base::failure(" NegfResult::saveTE(): Failed to open file " 
                        + fileName + ".");
            }

            // Energy
            out << "ENERGY" << endl;
            out << mE.n_elem << endl;
            out << mE;
            //k-points
            out << "KPOINTS" << endl;
            out << mk.n_rows << endl;
            if(!mk.is_empty()){
                out << mk;
            }
            // Transmission 
            if (mTE.isEnabled()){
                mTE.save(out, isText);
            }
            // Current
            for (int it = 0; it < mIop.size(); ++it){
                mIop[it].save(out, isText);
            }
            // Density of States
            if (mDOS.isEnabled()){
                mDOS.save(out, isText);
            }
            // Non-equilibrium electron density
            for (int it = 0; it < mnOp.size(); ++it){
                mnOp[it].save(out, isText);
            }        
            // Equilibrium electron density
            for (int it = 0; it < mpOp.size(); ++it){
                mpOp[it].save(out, isText);
            } 
        }else{
            
        }

    }
}

void CohRgfLoop::enableTE(uint N){
    mTE.tag = "TRANSMISSION";
    mTE.N = N;
}

void CohRgfLoop::enableI(uint N, uint ib, uint jb){
    mIop.push_back(RgfResult("CURRENT", N, ib, jb));
    mThisIop.push_back(cxmat_vec());
}

void CohRgfLoop::enableDOS(uint N){
    mDOS.tag = "DOS";
    mDOS.N = N;
}

void CohRgfLoop::enablen(uint N, int ib){    
    mnOp.push_back(RgfResult("n", N, ib, ib));
    mThisnOp.push_back(cxmat_vec());

}

void CohRgfLoop::enablep(uint N, int ib){    
    mpOp.push_back(RgfResult("p", N, ib, ib));
    mThispOp.push_back(cxmat_vec());

}

}
}

/**
 * Python exporters.
 */
namespace qmicad{
namespace python{
using namespace negf;

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enableTE, enableTE, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enableI, enableI, 0, 3)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enableDOS, enableDOS, 0, 1)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enablen, enablen, 0, 2)
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(CohRgfLoop_enablep, enablep, 0, 2)
void export_CohRgfLoop(){
    class_<CohRgfLoop, shared_ptr<CohRgfLoop> >("CohRgfLoop", 
            init<const Workers&, 
            optional<uint, double, dcmplx, bool, string> >())
        .def("run", &CohRgfLoop::run)
        .def("save", &CohRgfLoop::save)
        .def("enableTE", &CohRgfLoop::enableTE, CohRgfLoop_enableTE())
        .def("enableI", &CohRgfLoop::enableI, CohRgfLoop_enableI())
        .def("enableDOS", &CohRgfLoop::enableDOS, CohRgfLoop_enableDOS())
        .def("enablen", &CohRgfLoop::enablen, CohRgfLoop_enablen())
        .def("enablep", &CohRgfLoop::enablep, CohRgfLoop_enablep())
    ;
}

}
}



