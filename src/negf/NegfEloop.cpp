/* 
 * File:   NegfEloop.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 28, 2014, 7:20 PM
 */

#include "NegfEloop.h"

void NegfEloop::prepare() {
    if (mIAmMaster){
        cout << " NEGF: |"; 
    }
}


void NegfEloop::preCompute(int il){
    double E = mL(il);
    mnegf = shared_ptr<Negf>(new Negf(mnp, E));
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
void NegfEloop::computeTE(uword N){
    mTEproc.push_back(mnegf->TEop(N));
    //cout << "DBG: " << mTEl.back() << endl;
}

void NegfEloop::gather(vector<cx_mat> &Mproc, list<cx_mat> &M){
    if(!mIAmMaster){
        // slaves send their local data
        mpi::gather(mWorkers, Mproc, mMasterId);

    // The master collects data        
    }else{
        // Collect T(E)
        vector<vector<cx_mat> >Mncpu(mN);
        mpi::gather(mWorkers, Mproc, Mncpu, mMasterId);
        
        vector<vector<cx_mat> >::iterator it;
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

        list<cx_mat>::iterator it = mTE.begin();
        for (int i = 0; i < mL.N(); ++i, ++it){

            if (mTEn == 1){
                rowvec rvec(2); 
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



