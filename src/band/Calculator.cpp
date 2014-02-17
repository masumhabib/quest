/* 
 * File:   Calculator.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 * 
 * Created on April 23, 2013, 11:21 AM
 * 
 * Description: Base of all calculator classes.
 * 
 */

#include "Calculator.h"


const string Calculator::mOutExt = "dat";

#ifdef USE_MPI
Calculator::Calculator(mpi::communicator *world, Atoms* atm, TightBinder *tb, string prefix):
mWorld(world),
#else
Calculator::Calculator(Atoms* atm, TightBinder *tb, string prefix):
#endif
Printable(" " + prefix),
mSuperCell(atm), 
mTightBinder(tb)
{

    assert(mSuperCell);
    assert(mTightBinder);
    
#ifdef USE_MPI
    assert(mWorld);
    mMyCpuId = mWorld->rank();
    mnCpu = mWorld->size();
    mIamMaster = (mMasterId == mMyCpuId);
#else
    mMyCpuId = 0;
    mnCpu = 1;
    mIamMaster = true;
#endif
    
    mV = 0;
    mIsSaveNeigh = false;
    mTitle = "Calculator";

}

Calculator::~Calculator() {
}

void Calculator::generateNeighHam(cube &H, vector<svec>& r, 
        const Atoms* superCellI, const Atoms* superCellJ){

    // by default we use our own supercell
    if (superCellI == NULL){
        superCellI = mSuperCell;
    }
    if (superCellJ == NULL){
        superCellJ = superCellI;
    }
    
    // set up the nearest neighbors
    vector<lcoord> nc;      // Nearest neighbor coordinates
    int N = 2;              // Consider N nearest neighbors
    for (int i = -N; i <= N; ++i){
        for (int j = -N; j <= N; ++j){
            nc.push_back(lcoord(i,j,0));
        }
    }

    
    // CHECKS
    if (superCellI->NumOfAtoms() == 0 || superCellJ->NumOfAtoms() == 0){
        throw runtime_error(" No atoms found in the supercell.");
    }
    
    int noi = superCellI->NumOfOrbitals();     // number of orbitals
    int noj = superCellJ->NumOfOrbitals();     // number of orbitals
    lvec a = superCellI->LatticeVector();      // lattice vector

    //----- Construct the hamiltonian -----
    // loop through the nearest neighbors and get the matrix elements
    Atoms atmj;     // neighboring cell
    Atoms neighs;   // all neighbors
    H.set_size(noi, noj, nc.size());       // Hamiltonian matrix for the nearest neighbors
    for(int in = 0; in != nc.size(); ++in){
        // position vector of the neighbors
        r.push_back(nc[in].n1*a.a1 + nc[in].n2*a.a2 + nc[in].n3*a.a3);
        // generate nearest neighbors
        atmj = *superCellJ + nc[in];
        
        //generate the hamiltonian matrices
        H.slice(in) = mTightBinder->generateHamiltonian(*superCellI, atmj);

        // Do we want to generate neighboring cells?
        if(mIsSaveNeigh){
            neighs += atmj;            
            /*//DBG
            if (in == 0){
                neighs += atmj;
            } */  
            }  
        }
    
    if(mIsSaveNeigh && mIamMaster){
        // create directory if not exist
        int status = mkdir(mOutputPath.c_str(), S_IRWXU | S_IXGRP);
        if(status != 0 && errno != EEXIST){
            throw runtime_error("Cannot create/access directory " + mOutputPath);
        }
        
        string fileName = mOutputPath + "/" + mNeighFileName;
        neighs.exportGjf(fileName);
    }
}

void Calculator::generateHk(cx_mat& Hk, cube& H, rowvec& k, vector<svec>& r){
    //---- calculate H(k) ----
    // set Hk = H00 + i*0
    Hk.zeros();
    // OLD Hk.set_real(H.slice(0));
    dcmplx ith;
    
    //OLD for(int in = 1; in != H.n_slices; ++in){
    for(int in = 0; in != H.n_slices; ++in){

        ith = i*(r[in].col(X)*k(X)
              + r[in].col(Y)*k(Y));
        //OLD Hk = Hk + H.slice(in)*exp(ith) + trans(H.slice(in))*exp(-ith);
        Hk = Hk + H.slice(in)*exp(ith);
    }

}

void Calculator::readOptions(const TwOpts& opts){
    TwOpts::CalculationOpts cpts = opts.CalculationOptions();
    mOutputPath = cpts.OutPath;
    mNeighFileName = cpts.neigh.FileName;
    mIsSaveNeigh = (cpts.neigh.Generate == "true");
    mOutputFileName = cpts.FileName;        
}


string Calculator::toString() const { 
    stringstream ss;
    ss << mPrefix << " outputPath: " << mOutputPath << endl;
    ss << mPrefix << " genNeigh: " << mIsSaveNeigh << endl;
    ss << mPrefix << " neighFileName: " << mNeighFileName << endl;

    return ss.str(); 
};

void Calculator::assignCPUs(int N, int& myStartPoint, int& myEndPoint){
    int quo = N/mnCpu;
    int rem = N%mnCpu;
    myStartPoint = mMyCpuId*quo + (mMyCpuId < rem ? mMyCpuId:rem);
    myEndPoint = myStartPoint + quo + (mMyCpuId < rem ? 1:0) - 1;
}


void Calculator::generatePotential(const string& bias, const vector<TwOpts::Potential>& potential){
    
    mV = stod(bias);
    
    int na = mSuperCell->NumOfAtoms();
    mPotential.set_size(na);
    mPotential.fill(mV); // By default, all the atoms have same potential as bias
    
    for (int i = 0; i < potential.size(); ++i){
        
        
        int start = stoi(potential[i].StartAtom);
        int end = stoi(potential[i].EndAtom);
        double ratio = stod(potential[i].Ratio);
        
        if(start < 1 || end > na || start > end){ //check the range
            throw runtime_error(" Atom indices in <Potential> section is out of bounds.");  
        }
        mPotential.rows(start-1,end-1).fill(-mV*ratio);
    }
}


void Calculator::genKPoints(const string& kp, mat& k){
    stringstream ssk(trim(kp));                                                                                   

    // Read the k-points                            
    rowvec newk(2);
    // read new point
    ssk >> newk(X) >> newk(Y);                                                                
    // insert new point
    k.insert_rows(k.n_rows,newk);                                

}

void Calculator::genKLines(const TwOpts::KLine& kl, mat& k){ 
    
    // Get the starting ending and num of new k-points
    double kxs = stod(kl.KxStart);
    double kys = stod(kl.KyStart);
    double kxe = stod(kl.KxEnd);
    double kye = stod(kl.KyEnd);
    double delk = stod(kl.DelK);
    
    uword nknew = 1;
    if (abs(delk) > 0){
        double dk = sqrt(pow(kxs - kxe, 2) + pow(kys - kye,2));
        nknew = uword(dk/delk) + 1;
    }
    
    // Generate k points along the line                             
    mat newk(nknew,2);
    // generate kx and ky
    newk.col(X) = linspace<colvec>(kxs, kxe, nknew);
    newk.col(Y) = linspace<colvec>(kys, kye, nknew);
    // insert new k-points
    k.insert_rows(k.n_rows,newk);

}

void Calculator::genKSurface(const TwOpts::KSurface& ks, mat& k){                                                                                 

    double kxmin = stod(ks.KxStart);
    double kxmax = stod(ks.KxEnd);
    double kymin = stod(ks.KyStart);
    double kymax = stod(ks.KyEnd);
    double delk = stod(ks.DelK);
    
    uword nkx = 1;
    uword nky = 1;
    
    if (abs(delk) > 0){
        nkx = uword(abs(kxmax - kxmin)/delk) + 1;
        nky = uword(abs(kymax - kymin)/delk) + 1;
    }

    // Generate k mesh grid    
    rowvec newkx(nkx);
    colvec newky(nky);
    // generate kx and ky
    newkx = linspace<rowvec>(kxmin, kxmax, nkx);
    newky = linspace<colvec>(kymin, kymax, nky);
    // insert the new k-points
    for(int ikx = 0; ikx < nkx; ++ikx){
        mat newk(nky,2);

        newk.col(X).fill(newkx(ikx));
        newk.col(Y) = newky;
        k.insert_rows(k.n_rows,newk);
    }
}

void Calculator::genKParallelo(const TwOpts::KParallelo& kpar, mat& k){                                                                                 

    double k1x = stod(kpar.K1x);
    double k1y = stod(kpar.K1y);
    double k2x = stod(kpar.K2x);
    double k2y = stod(kpar.K2y);
    double delk = stod(kpar.DelK);
    
    double k1 = sqrt(k1x*k1x + k1y*k1y);
    double k2 = sqrt(k2x*k2x + k2y*k2y);
    
    uword nk1 = 1;
    uword nk2 = 1;
    
    if (abs(delk) > 0){
        nk1 = uword(k1/delk) + 1;
        nk2 = uword(k2/delk) + 1;
    }

    // Generate k mesh grid in a parallelogram    
    colvec newk1x(nk1);
    colvec newk1y(nk1);
    // generate kx and ky along vector k1
    newk1x = linspace<colvec>(0, k1x, nk1);
    newk1y = linspace<colvec>(0, k1y, nk1);
    // del k vector along k2
    k2x = k2x/k2*delk;
    k2y = k2y/k2*delk;
    // insert the new k-points
    for(int ik2 = 0; ik2 < nk2; ++ik2){
        mat newk(nk1,2);
        newk.col(X) = newk1x + k2x*ik2;
        newk.col(Y) = newk1y + k2y*ik2;
        k.insert_rows(k.n_rows,newk);
    }
}

void Calculator::genKSmartCircle(const TwOpts::KSmartCircle& kc, mat& k){
    
    double centerX = stod(kc.CenterX);
    double centerY = stod(kc.CenterY);
    double radiusStart = stod(kc.RadiusStart);
    double radiusEnd = stod(kc.RadiusEnd);
    double thetaStart = stod(kc.ThetaStart)*datum::pi/180.0;
    double thetaEnd = stod(kc.ThetaEnd)*datum::pi/180.0;
    double delk = stod(kc.DelK);
    
    radiusStart = abs(radiusStart);
    radiusEnd = abs(radiusEnd); 
    
    uword nRadius = 1;
    if (abs(delk) > 0){
        nRadius = uword(abs(radiusEnd - radiusStart)/delk) + 1;
    }
    
    colvec radiusGrid = linspace<colvec>(radiusStart, radiusEnd, nRadius);
    
    for(int it = 0; it < nRadius; ++it){
        uword thetaN = 1;
        if (abs(delk) > 0){
            thetaN = uword(abs(thetaEnd-thetaStart)*radiusGrid[it]/delk) + 1;
        }
        colvec thetaGrid = linspace<colvec>(thetaStart, thetaEnd, thetaN);
        
        mat newk;
        newk.set_size(thetaN, 2);
        newk.col(X) = radiusGrid(it)*cos(thetaGrid) + centerX;
        newk.col(Y) = radiusGrid(it)*sin(thetaGrid) + centerY;
        
        k.insert_rows(k.n_rows, newk);
    }
}

void Calculator::genKCircle(const TwOpts::KCircle& kc, mat& k){
    
    double centerX = stod(kc.CenterX);
    double centerY = stod(kc.CenterY);
    double radiusStart = stod(kc.RadiusStart);
    double radiusEnd = stod(kc.RadiusEnd);
    double delRadius = stod(kc.DelRadius);
    double thetaStart = stod(kc.ThetaStart)*datum::pi/180.0;
    double thetaEnd = stod(kc.ThetaEnd)*datum::pi/180.0;
    double delTheta = stod(kc.DelTheta)*datum::pi/180.0;
    
    radiusStart = abs(radiusStart);
    radiusEnd = abs(radiusEnd); 
    uword nRadius = 1;
    if (abs(delRadius) > 0){
        nRadius = uword(abs(radiusEnd - radiusStart)/delRadius) + 1;
    }    
    colvec radiusGrid = linspace<colvec>(radiusStart, radiusEnd, nRadius);

    uword thetaN = 1;    
    if (abs(delTheta) > 0){
        thetaN = uword(abs(thetaEnd-thetaStart)/delTheta) + 1;
    }
    colvec thetaGrid = linspace<colvec>(thetaStart, thetaEnd, thetaN);
    
    for(int it = 0; it < nRadius; ++it){
        
        mat newk;
        newk.set_size(thetaN, 2);
        newk.col(X) = radiusGrid(it)*cos(thetaGrid) + centerX;
        newk.col(Y) = radiusGrid(it)*sin(thetaGrid) + centerY;
        
        k.insert_rows(k.n_rows, newk);
    }
}


string Calculator::Bias() const{
    stringstream ssbias;
    ssbias.precision(4);
    ssbias << std::fixed << mV;
    
    return ssbias.str();
}
