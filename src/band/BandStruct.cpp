/* 
 * File:   BandStruct.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 6, 2013, 9:25 AM
 */

#include "band/BandStruct.h"

namespace qmicad{
namespace band{

BandStruct::BandStruct(shared_ptr<mat> pk, const BandStructParams &bp, 
        const Workers &workers, bool saveAscii): ParLoop(workers, pk->n_rows), 
        mp(bp), mk(pk) , mSaveAscii(saveAscii), mCalcEigV(false),
        mbar("  EK: ",  pk->n_rows)
{    
    
    mlb = mp.ne/2 - mp.nb/2;     // lowest band to calculate
    mub = mp.ne/2 + mp.nb/2-1;   // highest band to calculate  

    if (mub < mlb){
        swap(mub, mlb);
    }
    
    if (mlb < 0L){
        mlb = 0;
    }
    
    if (mub >= mp.no){
        mub = mp.no-1;
    }    
}

void BandStruct::prepare(){
    
    //Runtime checks.
    if (mp.no != mp.H(0)->n_rows){
        stringstream out;
        out << "Size of Hamiltonian matrix does not match with number "
            << "of orbitals: no = " << mp.no << " size(H) = " << mp.H(0)->n_rows
            << "x" << mp.H(0)->n_cols << "";
        throw runtime_error(out.str());
    }
        
    // Setup
    mThisE.set_size(mMyN, mub - mlb + 1);     // Eigen energy: (# of kpoints)  x  (# of bands + size of k vector).
    
    mWorkers.Comm().barrier();
    mbar.start();
}

void BandStruct::preCompute(int il){
    
}

void BandStruct::postCompute(int il){
    ++mbar;
}

void BandStruct::compute(int il){
    
    uint    N = mp.no;
    cxmat   Hk(N, N, fill::zeros);
    row k = mk->row(il);
    // loop over half the nearest neighbors
    for (uint ih = 0; ih < mp.H.n_elem; ++ih){
        lcoord lc = mp.lc(ih);
        if (lc.n1 == 0 && lc.n2 == 0 && lc.n3 ==0){
            Hk += *mp.H(ih);
        }else{
            double th = dot(k, mp.lv*lc);    // th = k.*(R*n)
            Hk += (*mp.H(ih))*exp(i*th) + trans(*mp.H(ih))*exp(-i*th);
        } 
    }
    
    if (mCalcEigV){
        //Eigen vector calculation.
    }else{
        col E = eig_sym(Hk);
        int ik = il - mMyStart;
        mThisE.row(ik) = trans(sort(E.rows(mlb, mub)));
    }
}

void BandStruct::collect(){
    mWorkers.Comm().barrier();
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
        // open file.
        ofstream out;
        if(mSaveAscii){
            out.open(fileName.c_str(), ostream::binary);
        }else{
            out.open(fileName.c_str());
        }
        if (!out.is_open()){
            throw ios_base::failure(" NegfResult::saveTE(): Failed to open file " 
                    + fileName + ".");
        }
        
        out << "EK" << endl;    // tag
        out << mk->n_rows << endl; // # of k points
        out << 1 << " " << mE.n_cols << endl; // 1xnb matrix; nb = number of bands

        for (int ik = 0; ik < mk->n_rows; ++ik){
            out << mk->row(ik);
            out << mE.row(ik);
        }

        out.close();
        
        if (mCalcEigV){
        }

    }
}


}
}

