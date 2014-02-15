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
        mTE("TE", 1, saveAscii),
        mI1op("I1", 0, saveAscii), mI1N(0), mI1sxN(0), mI1syN(0), mI1szN(0)
{
    mprog   = 0;
    mprogmx = 80;
}
    

void NegfEloop::prepare() {
    if (mIAmMaster){
        cout << "  NEGF: |"; 
    }
}

void NegfEloop::preCompute(int il){
    double E = mL(il);
    mnegf = shared_ptr<CohRgfa>(new CohRgfa(mnp, E));
};

void NegfEloop::compute(int il){
    negfresult r;
    r.first = mL(il);               // this energy
    // Transmission
    if(mTE.isEnabled()){
        r.second = mnegf->TEop(mTE.N);  // Calculate  T(E)
        mThisTE.push_back(r);  
    }
    // Current
    if(mI1op.isEnabled()){
        r.second = mnegf->I1Op(mI1op.N);
        mThisI1op.push_back(r);
    }
}

void NegfEloop::postCompute(int il){
    mnegf.reset();          // free up memory
    stepCompleted();        // Show feedback
}

void NegfEloop::collect(){
    
    // Transmission
    if(mTE.isEnabled()){
        gather(mThisTE, mTE);
    }
    
    // Current
    if(mI1op.isEnabled()){
        gather(mThisI1op, mI1op);
    }
    
    if (mIAmMaster){
        cout << "|"; 
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
        // Transmission 
        if (mTE.isEnabled()){
            mTE.save(fileName);
        }
        
        // Current
        if (mI1op.isEnabled()){ 
            NegfResultList  I1("I1", mI1N, mI1op.saveAscii);
            NegfResultList  I1sx("I1sx", mI1sxN, mI1op.saveAscii);
            NegfResultList  I1sy("I1sy", mI1syN, mI1op.saveAscii);
            NegfResultList  I1sz("I1sz", mI1szN, mI1op.saveAscii);
            
            list<negfresult> &R = mI1op.R;
            list<negfresult>::iterator it;
            for (it = R.begin(); it != R.end(); ++it){
                negfresult r;
                r.first = it->first;
                // Charge current
                if(mI1N){
                    r.second = trace(it->second, mI1N);
                    I1.R.push_back(r);
                }
                
                // Spin currents
                if(mI1sxN){
                    r.second = trace(trace(it->second, 2)*sx());
                    I1sx.R.push_back(r);
                }
                if(mI1syN){
                    r.second = trace(trace(it->second, 2)*sy());
                    I1sy.R.push_back(r);
                }
                if(mI1szN){
                    r.second = trace(trace(it->second, 2)*sz());
                    I1sz.R.push_back(r);
                }
            }
            if(mI1N){
                I1.save(fileName);
            }
            if(mI1sxN){
                I1sx.save(fileName);
            }
            if(mI1syN){
                I1sy.save(fileName);
            }
            if(mI1szN){
                I1sz.save(fileName);
            }
        }
    }
}

void NegfEloop::enableTE(uint N){
    mTE.N = N;
}

void NegfEloop::disableTE(){
    mTE.N = 0;
}

void NegfEloop::enableI1(uint N){
    mI1N = N;
    if (mI1N > mI1op.N){
        mI1op.N = mI1N;
    }
}

void NegfEloop::enableI1sx(uint N){
    mI1sxN = N;
    if (mI1sxN > 0 && mI1op.N < 2){
        mI1op.N = 2;
    }    
}

void NegfEloop::enableI1sy(uint N){
    mI1syN = N;
    if (mI1syN > 0 && mI1op.N < 2){
        mI1op.N = 2;
    }    
}

void NegfEloop::enableI1sz(uint N){
    mI1szN = N;
    if (mI1szN > 0 && mI1op.N < 2){
        mI1op.N = 2;
    }    
}


}
