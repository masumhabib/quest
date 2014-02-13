/* 
 * File:   NegfEloop.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "NegfEloop.h"
#include "NegfResult.h"

namespace qmicad{

NegfEloop::NegfEloop(const VecGrid &E, const NegfParams &np, 
        const mpi::communicator &workers, bool saveAscii):
        ParLoop<double>(E, workers), mnp(np), 
        mTE("TE", true, saveAscii)
{
    mprog   = 0;
    mprogmx = 80;
}
    

void NegfEloop::prepare() {
    if (mIAmMaster){
        cout << " NEGF: |"; 
    }
}


void NegfEloop::preCompute(int il){
    double E = mL(il);
    mnegf = shared_ptr<CohRgfa>(new CohRgfa(mnp, E));
};

void NegfEloop::compute(int il){
    if(mTE.enabled){
        negfresult r;
        r.second = mnegf->TEop(mTE.N);  // Calculate  T(E)
        r.first = mL(il);               // this energy
        mThisTE.push_back(r);  
    }
}

void NegfEloop::postCompute(int il){
    mnegf.reset();          // free up memory
    stepCompleted();        // Show feedback
}

void NegfEloop::collect(){
    
    if(mTE.enabled){
        gather(mThisTE, mTE);   // Gather T(E)
    }
    
    if (mIAmMaster){
        cout << "|" << endl; 
    }
}

void NegfEloop::stepCompleted(){
    if (++mprog * mprogmx/mN != 0){
        cout << "*";
        mprog = 0;
    }
}

void NegfEloop::gather(vector<negfresult> &thisR, NegfResultList &all){
    if(!mIAmMaster){
        // slaves send their local data
        mpi::gather(mWorkers, thisR, mMasterId);

    // The master collects data        
    }else{
        // Collect T(E)
        vector<vector<negfresult> >gatheredR(mN);
        mpi::gather(mWorkers, thisR, gatheredR, mMasterId);
        
        // merge and store results on mTE list.
        vector<vector<negfresult> >::iterator it;
        for (it = gatheredR.begin(); it != gatheredR.end(); ++it){
            all.R.insert(all.R.end(), it->begin(), it->end());
        }
        // Sort the results based on energy.
        all.sort();
    }       
}

void NegfEloop::save(string fileName){
    if(mIAmMaster){
        if (mTE.enabled){
            mTE.save(fileName);
        }
    }
}

void NegfEloop::enableTE(uint N){
    mTE.enabled = true;
    mTE.N = N;
}

void NegfEloop::disableTE(){
    mTE.enabled = false;
}


}
