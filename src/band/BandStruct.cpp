/* 
 * File:   BandStruct.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 * 
 * Created on April 6, 2013, 9:25 AM
 */

#include "BandStruct.h"

namespace qmicad{

BandStruct::BandStruct(shared_ptr<mat> pk, const BandStructParams &bp, 
        const Workers &workers, bool saveAscii): ParLoop(workers, pk->n_rows), 
        mSaveAscii(saveAscii), mp(bp), mk(pk) , mbar("  NEGF: ",  pk->n_rows),
        mCalcEigV(false)
{    
    mlb = mp.ne/2 - mp.nb/2;     // lowest band to calculate
    mub = mp.ne/2 + mp.nb/2-1;   // highest band to calculate    
}

void BandStruct::prepare(){
    
    mThisE.set_size(mMyN, mk->n_cols + mp.nb);     // Eigen energy: (# of kpoints)  x  (# of bands + size of k vector).
    mbar.start();
}

void BandStruct::preCompute(int il){
    
}

void BandStruct::postCompute(int il){
    ++mbar;
}

void BandStruct::compute(int il){
    
    uint    N = mp.H(0)->n_rows;
    cxmat   Hk(N, N, fill::zeros);
    row k = mk->row(il);
    // loop over half the nearest neighbors
    for (int ih = 0; ih < mp.H.n_elem; ++ih){
        lcoord lc = mp.lc(ih);

        if (lc.n1 == 0 && lc.n2 == 0 && lc.n3 ==0){
            Hk += *mp.H(ih);
        }else{
            double th = dot(k, mp.lv*lc);    // th = k.*(R*n)
            Hk += (*mp.H(ih))*exp(th) + trans(*mp.H(ih))*exp(-th);
        }
    }
    if (mCalcEigV){
        //Eigen vector calculation.
    }else{
        col E = eig_sym(Hk);
        int ik = il - mMyStart;
        mThisE.row(ik) = join_rows(k, sort(trans(E.rows(mlb, mub))));
    }
}

void BandStruct::collect(){
    mbar.complete();
    
    if (mCalcEigV){
        
    }else{
        // Gather data from all the processes.
        if(!mWorkers.IAmMaster()){    
            // slaves send their local data
            mpi::gather(mWorkers.Comm(), mThisE, mWorkers.MasterId());

        // The master collects data        
        }else{
            // Collect T(E)
            vector<mat> gatheredE(mN);
            mpi::gather(mWorkers.Comm(), mThisE, gatheredE, mWorkers.MasterId());

            // merge and store results on mTE list.
            vector<mat>::iterator it;
            for (it = gatheredE.begin(); it != gatheredE.end(); ++it){
                mE.insert_rows(mE.n_rows, *it);
            }
        }
    }
}

void BandStruct::save(string fileName){
    if(mWorkers.IAmMaster()){
        if (mCalcEigV){
        }else{
            fileName = fileName + "EK" + ".dat";
            
            // open file.
            ofstream outFile;
            if(mSaveAscii){
                outFile.open(fileName.c_str(), ostream::binary);
            }else{
                outFile.open(fileName.c_str());
            }
            if (!outFile.is_open()){
                throw ios_base::failure(" NegfResult::saveTE(): Failed to open file " 
                        + fileName + ".");
            }
            
            // save to the file
            outFile << mE;
            outFile.close();
        }
    }
}


}