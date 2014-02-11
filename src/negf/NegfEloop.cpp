/* 
 * File:   NegfEloop.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "NegfEloop.h"

namespace qmicad{
    

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
    if(mCalcType & TE){
        computeTE(mTEn);
    }
}

void NegfEloop::postCompute(int il){
    mnegf.reset(); // free up memory
    stepCompleted();
}

void NegfEloop::collect(){
    
    if(mCalcType & TE){
        collectTE();
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
void NegfEloop::computeTE(uint N){
    mTEproc.push_back(mnegf->TEop(N));
    //cout << "DBG: " << mTEl.back() << endl;
}

void NegfEloop::gather(vector<cxmat> &Mproc, list<cxmat> &M){
    if(!mIAmMaster){
        // slaves send their local data
        mpi::gather(mWorkers, Mproc, mMasterId);

    // The master collects data        
    }else{
        // Collect T(E)
        vector<vector<cxmat> >Mncpu(mN);
        mpi::gather(mWorkers, Mproc, Mncpu, mMasterId);
        
        vector<vector<cxmat> >::iterator it;
        for (it = Mncpu.begin(); it != Mncpu.end(); ++it){
            M.insert(M.end(), it->begin(), it->end());
        }
    }       
}

void NegfEloop::collectTE(){
    gather(mTEproc, mTE);
}

void NegfEloop::saveTE(string fileName){
    if(mIAmMaster){
        // save to a file
        ofstream transFile(fileName.c_str());
        if (!transFile.is_open()){
            throw ios_base::failure(" NegfLoop::saveTE(): Failed to open file " 
                    + fileName + ".");
        }

        list<cxmat>::iterator it = mTE.begin();
        for (int i = 0; i < mL.N(); ++i, ++it){

            if (mTEn == 1){
                row rvec(2); 
                rvec(0) = mL(i);
                rvec(1) = real((*it)(0,0));
                transFile << rvec;  
            }else{
                transFile << "E = " << mL(i) << "TE: " << endl << *it << endl;
            }
        }

        transFile.close();
    }    
}

}

