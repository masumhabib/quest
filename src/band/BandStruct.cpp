/* 
 * File:   BandStruct.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 6, 2013, 9:25 AM
 */

#include "band/BandStruct.h"

namespace qmicad{
namespace band{

BandStruct::BandStruct(const Workers &workers, uint nn,  bool orthoBasis, 
        bool calcEigV, const string &prefix): Printable(prefix), 
        mWorkers(workers), mnn(nn),  mOrthoBasis(orthoBasis), mCalcEigV(calcEigV), 
        mbar("  EK: "), mH(nn), mS(nn), mlc(nn)
{    
    mTitle = "Band Structure";
}

void BandStruct::nb(uint nb){
    mnb = nb;
}

void BandStruct::ne(uint ne){
    mne = ne;
}

void BandStruct::lv(const lvec& lv){
    mlv = lv;
}

void BandStruct::lc(field<lcoord> lc){
    if (lc.n_elem != mnn){
        throw runtime_error("In CohRgfLoop::lc() size of lc should be equal to number of neighbors.");
    }    
    
    mlc = lc;
}

void BandStruct::lc(const lcoord &lc, int nn){
    mlc(nn) = lc;
}

void BandStruct::k(const mat& k){
    mk = k;
    mN = mk.n_rows;
    mWorkers.assignCpus(mMyStart, mMyEnd, mMyN, mN);
    mbar.expectedCount(mN);
}

void BandStruct::H(const field<shared_ptr<cxmat> >& H)
{
    if (H.n_rows != mnn){
        throw runtime_error("In CohRgfLoop::H() size of H does not match with number of neighbors");
    }
    
    mH = H;
}

void BandStruct::S(const field<shared_ptr<cxmat> >& S)
{
    if (S.n_rows != mnn){
        throw runtime_error("In CohRgfLoop::S() size of S does not match with number of neighbors");
    }
    
    mS = S;
}


void BandStruct::H(shared_ptr<cxmat> H, int ineigh){
    mH(ineigh) = H;
    mno = H->n_rows;
}

void BandStruct::S(shared_ptr<cxmat> S, int ineigh){
    mS(ineigh) = S;
}


void BandStruct::run(){
    prepare();
    for(int il = mMyStart; il <= mMyEnd; ++il){
        preCompute(il);
        compute(il);
        postCompute(il);
    }    
    collect();
}

void BandStruct::save(string fileName, bool saveAsText){
    if(mWorkers.IAmMaster()){
        if (saveAsText){ 
            // open file.
            ofstream out;
            out.open(fileName.c_str(), ios::app);
            
            if (!out.is_open()){
                throw ios_base::failure(" NegfResult::saveTE(): Failed to open file " 
                        + fileName + ".");
            }

            out << "EK" << endl;    // tag
            out << mk.n_rows << endl; // # of k points
            out << 1 << " " << mE.n_cols << endl; // 1xnb matrix; nb = number of bands

            for (int ik = 0; ik < mk.n_rows; ++ik){
                out << mk.row(ik);
                out << mE.row(ik);
            }

            out.close();

            if (mCalcEigV){
            }
        }else{ // binary file
            
        }
    }
}

void BandStruct::resetBandIndices(){
    if (mnb > mno || mnb < 0){
        mnb = mno;
    }
    
    mlb = mne/2 - mnb/2;     // lowest band to calculate
    mub = mne/2 + mnb/2-1;   // highest band to calculate  

    if (mub < mlb){
        swap(mub, mlb);
    }
    
    if (mlb < 0L){
        mlb = 0;
    }
    
    if (mub >= mno){
        mub = mno-1;
    } 
}

void BandStruct::prepare(){
    resetBandIndices();
    //Runtime checks.
    if (mno != mH(0)->n_rows){
        stringstream out;
        out << "Size of Hamiltonian matrix does not match with number "
            << "of orbitals: no = " << mno << " size(H) = " << mH(0)->n_rows
            << "x" << mH(0)->n_cols << "";
        throw runtime_error(out.str());
    }
        
    // Setup
    mThisE.set_size(mMyN, mub - mlb + 1);     // Eigen energy: (# of kpoints)  x  (# of bands + size of k vector).
    
    mWorkers.Comm().barrier();
    mbar.start();
}

void BandStruct::preCompute(long il){
    
}

void BandStruct::postCompute(long il){
    ++mbar;
}

void BandStruct::compute(long il){
    
    uint    N = mno;
    cxmat   Hk(N, N, fill::zeros);
    row k = mk.row(il);
    // loop over the nearest neighbors
    for (uint ih = 0; ih < mH.n_elem; ++ih){
        lcoord lc = mlc(ih);
        if (lc.n1 == 0 && lc.n2 == 0 && lc.n3 ==0){
            Hk += *mH(ih);
        }else{
            double th = dot(k, mlv*lc);    // th = k.*(R*n)
            Hk += (*mH(ih))*exp(i*th) + trans(*mH(ih))*exp(-i*th);
        } 
    }
    
    if (mCalcEigV){
        //Eigen vector calculation.
    }else{
        col E = eig_sym(Hk);
        long ik = il - mMyStart;
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




}}

