/* 
 * File:   BandStruct.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 * 
 * Created on April 6, 2013, 9:25 AM
 */

#include "BandStruct.h"

namespace qmicad{

BandStruct::BandStruct(shared_ptr<mat> pk, const BandStructParams &bp, 
        const Workers &workers): ParLoop(workers, pk->n_rows), 
        mp(bp), mk(pk)
{    
    mlb = mp.Ne/2 - mp.Nb/2;     // lowest band to calculate
    mub = mp.Ne/2 + mp.Nb/2-1;   // highest band to calculate    
}

void BandStruct::prepare(){
    
    mThisE.set_size(mMyN, mp.Nb);     // Eigen energy: (# of kpoints)  x  (# of bands).
    mbar.start();
}

void BandStruct::preCompute(int il){
    
}

void BandStruct::postCompute(int il){
    
}

void BandStruct::compute(int il){
    
    uint    N = mp.H(0)->n_rows;
    cxmat   Hk(N, N, fill::zeros);
    
    mat &k = *mk;
    double kx = k(il, X);
    double ky = k(il, Y);
    //double kz = k(il, Z);
    
    for (int ih = 0; ih < mp.H.n_elem; ++ih){
        double th = kx*mp.pv(ih)(X) 
                  + ky*mp.pv(ih)(Y);
                  //+ kz*mp.pv(ih)(Z);
        
        Hk += *mp.H(ih)*exp(th);
    }
    
    int ik = il - mMyStart;
    col E = eig_sym(Hk);
    mThisE.row(ik) = trans(E.rows(mlb, mub));
    
    dout << dbg << endl << mThisE;
//    // CHECKS
//    if(mk.empty()){
//        throw runtime_error(" No k-points found.");
//    }
//    
//    int Nk = mk.n_rows; // total number of k-points
//    ConsoleProgressBar progress(Nk,0," EK: ");
//    
//    // if number of bands is out of range then save all the bands
//    if (mNumBands <= 0 || mNumBands > mSuperCell->NumOfOrbitals()){
//        mNumBands = mSuperCell->NumOfOrbitals();
//    }
//        
//// show progress
//    if (mIamMaster){ 
//        progress.start();
//        ++progress;
//    }
//    
//    //----- Construct the hamiltonian -----
//    cube H;                     // Hamiltonian matrix for the nearest neighbors
//    cx_mat Hk(mSuperCell->NumOfOrbitals(), mSuperCell->NumOfOrbitals());       // H(k)
//    vector<svec> r;             // position vectors of the neighbors
//    generateNeighHam(H, r);
//    
//    // show progress
//    if (mIamMaster){ 
//        ++progress;
//    }
//    
//    //----- Solve for eigen energy ------
//    // loop through the k-points
//    //k-points assignment to CPUs
//    int myKStarts, myKEnds;
//    assignCPUs(Nk, myKStarts, myKEnds);
//
//#ifdef USE_MPI
//    //For showing nice progress bar that includes progress from
//    //all CPUs.
//    int myProgress = 0;
//    int yourProgress = 0;
//    bool recvInProgress = false;
//    mpi::request recvReq;
//#endif
//    
//    //Apply the potential
//    H.slice(0).diag() += mPotential;
//   
//    //---- Main calculation ----
//    int lb = mSuperCell->NumOfElectrons()/2 - mNumBands/2;     // lowest band to calculate
//    int ub = mSuperCell->NumOfElectrons()/2 + mNumBands/2-1;   // highest band to calculate    
//    mE.set_size(ub-lb+1, Nk);
//    mE.zeros();
//    for(int ik = myKStarts; ik <= myKEnds; ++ik){
//        //---- calculate H(k) ----
//        rowvec ki = mk(ik, span::all);
//        generateHk(Hk, H, ki, r);
//
//        // Solve eigenvalue
//        colvec E = eig_sym(Hk);
//        mE.col(ik) = E.rows(lb, ub);
//        
//#ifdef USE_MPI        
//        // The master collects progress report
//        // and shows the porgress bar. To reduce MPI communication, we'll only
//        // send our progress report if progress is at least 1%.
//        if(mIamMaster){
//            // We dont want to call reccv if one recv is in progress.
//            if(recvInProgress){
//                // update progress report if we receive report from slaves. 
//                if(recvReq.test()){
//                    progress += yourProgress;
//                    recvInProgress = false;
//                }
//            }else{
//                recvReq = mWorld->irecv(mpi::any_source, PROGRESS, yourProgress);
//                recvInProgress = true;
//            }  
//            ++progress;
//        }else{
//            // If we have made significant progress then let the master know.
//            if (++myProgress*100/Nk/mnCpu > 0){ 
//                mWorld->isend(mMasterId, PROGRESS, myProgress);
//                myProgress = 0;
//            }
//        }
//#else
//        ++progress;
//#endif
//    }
//    
//    // Collect eigenvalues from all of the CPUs. 
//    if(mIamMaster){
//        progress.complete();
//#ifdef USE_MPI
//        cout << endl; 
//        cout << " EK: Collecting data ... ";
//        mpi::reduce(*mWorld, mE, mE, plus<mat>(), mMasterId);
//        cout << "done.";
//        
//    }else{
//        mpi::reduce(*mWorld, mE, plus<mat>(), mMasterId);
//#endif        
//    }    
}

void BandStruct::collect(){
    
}

bool BandStruct::save(){
//    mat ek = mk.t();
//    ek.insert_rows(2,mE);
//    
//    // create directory if not exist
//    int status = mkdir(mOutputPath.c_str(), S_IRWXU | S_IXGRP);
//    if(status != 0 && errno != EEXIST){
//        throw runtime_error("Cannot create/access directory " + mOutputPath);
//    }
//    
//    string fileName = mOutputPath + "/" + mOutputFileName 
//                      + "_" + Bias() + "." + mOutExt;
//    return ek.save(fileName, raw_ascii);
}


}