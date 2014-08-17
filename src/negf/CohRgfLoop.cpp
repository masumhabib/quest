/* 
 * File:   NegfEloop.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "negf/CohRgfLoop.h"
#include "negf/RgfResult.h"

namespace qmicad{
namespace negf{

CohRgfLoop::CohRgfLoop(const Workers &workers, uint nb, double kT, dcmplx ieta, 
        bool orthogonal, string newprefix): Printable(newprefix), 
        mrgf(nb, kT, ieta, orthogonal, " " + newprefix), mbar("  NEGF: "),
        mWorkers(workers)
{
    mH0.set_size(nb, 1);
    mS0.set_size(nb, 1);
    mHl.set_size(nb+1, 1);
    mSl.set_size(nb+1, 1);
    mV.set_size(nb);
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

void CohRgfLoop::numTransNeighbors(uint n){
    uint nb = mrgf.nb();
    mH0.set_size(nb, n);
    mS0.set_size(nb, n);
    mHl.set_size(nb+1, n);
    mSl.set_size(nb+1, n);
    
    mpv0.set_size(nb, n);
    mpvl.set_size(nb+1, n);
}

void CohRgfLoop::H(const field<shared_ptr<cxmat> >& H0, 
        const field<shared_ptr<cxmat> >& Hl)
{
    if (H0.n_rows != mrgf.nb() || Hl.n_rows != mrgf.nb()+1){
        throw runtime_error("In CohRgfLoop::H() size of H0 or Hl does not match with number of blocks");
    }
    
    mH0 = H0;
    mHl = Hl;
}

void CohRgfLoop::S(const field<shared_ptr<cxmat> >& S0, 
        const field<shared_ptr<cxmat> >& Sl)
{
    if (S0.n_rows != mrgf.nb()){
        throw runtime_error("In CohRgfLoop::S() size of S0 does not match with number of blocks");
    }
    
    mS0 = S0;
    mSl = Sl;
}

void CohRgfLoop::V(const field<shared_ptr<vec> >& V){
    if (V.n_rows != mrgf.nb()){
        throw runtime_error("In CohRgfLoop::V() size of V does not match with number of blocks");
    }
    mV = V;
}

void CohRgfLoop::pv(const field<vec>& pv0, const field<vec>& pvl){
    mpv0 = pv0;
    mpvl = pvl;
}

void CohRgfLoop::H0(shared_ptr<cxmat> H0, int ib, int ineigh){
    mH0(ib, ineigh) = H0;
}

void CohRgfLoop::S0(shared_ptr<cxmat> S0, int ib, int ineigh){
    mS0(ib, ineigh) = S0;
}

void CohRgfLoop::Hl(shared_ptr<cxmat> Hl, int ib, int ineigh){
    mHl(ib, ineigh) = Hl;
}

void CohRgfLoop::Sl(shared_ptr<cxmat> Sl, int ib, int ineigh){
    mSl(ib, ineigh) = Sl;
}

void CohRgfLoop::V(shared_ptr<vec> V, int ib){
    mV(ib) = V;
}

void CohRgfLoop::pv0(const vec& pv0, int ib, int ineigh){
    mpv0(ib, ineigh) = pv0;
}

void CohRgfLoop::pvl(const vec& pvl, int ib, int ineigh){
    mpvl(ib, ineigh) = pvl;
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

string CohRgfLoop::toString() const {
    stringstream out;
    out << mrgf;

    return out.str();
}

void CohRgfLoop::run(){
    
    prepare();

    long n = npoints();
    long nE = mE.n_rows;
    long nk = mk.n_rows;
    uint nb = mrgf.nb();
    
    long myStart, myEnd, myN;
    // Assign E and k points to CPUs 
    mWorkers.assignCpus(myStart, myEnd, myN, n);
    // Loop over problem assigned to this CPU.
    long ik, ikPrev = -1, iE;
    for(long it = myStart; it <= myEnd; ++it){
        
        ik = it/nE;
        iE = it%nE;

        if (ik != ikPrev){ // change H and S matrices only for new k vectors.
            // Hamiltonian and overlap matrices. 
            field<shared_ptr<cxmat> > H0(nb);
            field<shared_ptr<cxmat> > S0(nb);
            field<shared_ptr<cxmat> > Hl(nb+1);
            field<shared_ptr<cxmat> > Sl(nb+1);

            if (nk != 0){ // Do a k-loop
                //@TODO: Needs memory optimization.

                vec k = mk.row(ik);

                // Calculate block diagonal matrices.
                for(int ib = 0; ib < nb; ++ib){
                    shared_ptr<cxmat> H0k = make_shared<cxmat>(mH0(ib,0)->n_rows, mH0(ib,0)->n_cols, fill::zeros);
                    shared_ptr<cxmat> S0k = make_shared<cxmat>(mS0(ib,0)->n_rows, mS0(ib,0)->n_cols, fill::zeros);

                    double th;
                    dcmplx expith;
                    for(int in = 0; in < mH0.n_cols; ++in){
                        th = dot(k, mpv0(ib, in));
                        expith = exp(i*th);
                        (*H0k) = (*H0k) + (*mH0(ib, in))*expith;                    
                        (*S0k) = (*S0k) + (*mS0(ib, in))*expith;
                    }
                    H0(ib) = H0k;
                    S0(ib) = S0k;
                }
                // Calculate lower block diagonals
                for(int ib = 0; ib <= nb; ++ib){
                    shared_ptr<cxmat> Hlk = make_shared<cxmat>(mHl(ib,0)->n_rows, mHl(ib,0)->n_cols, fill::zeros);
                    shared_ptr<cxmat> Slk = make_shared<cxmat>(mSl(ib,0)->n_rows, mSl(ib,0)->n_cols, fill::zeros);

                    double th;
                    dcmplx expith;
                    for(int in = 0; in < mHl.n_cols; ++in){
                        th = dot(k, mpvl(in));
                        expith = exp(i*th);
                        (*Hlk) = (*Hlk) + (*mHl(ib, in))*expith;
                        if (!mrgf.OrthoBasis()){
                            (*Slk) = (*Slk) + (*mSl(ib, in))*expith;
                        }
                    }
                    Hl(ib) = Hlk;
                    if (!mrgf.OrthoBasis()){
                        Sl(ib) = Slk;
                    }
                }
            }else{ // Do only E loop
                H0 = mH0.col(0);
                S0 = mS0.col(0);
                Hl = mHl.col(0);
                Sl = mSl.col(0);   
            }
            // set H and S.
            mrgf.H(H0, Hl);
            mrgf.S(S0, Sl);
            mrgf.V(mV);
            
            ikPrev = ik;
        }
        
        // set E
        mrgf.E(mE[iE]);
        
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


}
}




