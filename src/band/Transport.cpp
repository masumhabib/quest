/* 
 * File:   Transport.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 * 
 * Created on April 22, 2013, 8:05 PM
 */

#include "Transport.h"

#ifdef USE_MPI
Transport::Transport(mpi::communicator *world, Atoms* atm, TightBinder *tb, string prefix):
Calculator(world, atm, tb, prefix),
#else
Transport::Transport(Atoms* atm, TightBinder *tb, string prefix):
Calculator(atm, tb, prefix),
#endif
mnKpoints(0)
{
    mTitle = "Transport";
    
}

Transport::~Transport() {
}

void Transport::computeTek(){
    if (mDeviceType == "NonRGF"){
        computeTekNonRGF();
    }else if(mDeviceType == "RGF"){
        
    }
}


void Transport::computeTekNonRGF(){
    int NE = mE.n_cols;
    ConsoleProgressBar progress(NE+1, 0, " TE: ");
         
    // show progress
    if (mIamMaster){ 
        progress.start();
    }
    
    // ----- Construct the Hamiltonian matrices -----
    
    cube HB;                     // Bottom contact
    cube HT;                     // Top contact
    cube TBD;                    // <Bottom|H|Device>
    cube TDT;                    // <Device|H|Top>        
    cube HD;                     // Device
    vector<svec> rd;             // position vectors of the neighbors
    
    genHam(HB, TBD, HD, TDT, HT, rd); // get the Hamiltonians
    
    HB.slice(0).diag() += mPotential.elem(mBottomContact); // apply potential
    HT.slice(0).diag() += mPotential.elem(mTopContact);    // apply potential
    HD.slice(0).diag() += mPotential.elem(mDevice); //Apply the potential
    
    // show progress
    if (mIamMaster){ 
        ++progress;
    }
    
    // ----- Energy loop ------
    mTek.set_size(NE);
    cx_mat HkD;     // H(k) for device
    cx_mat HkB;     // H(k) for the bottom contact
    cx_mat HkT;     // H(k) for the top contact
    cx_mat TkBD;    // H(k) for <Bottom|H|Device>
    cx_mat TkDT;    // H(k) for <Device|H|Top>

    HkD.copy_size(HD.slice(0));
    HkB.copy_size(HB.slice(0));
    HkT.copy_size(HT.slice(0));
    TkBD.copy_size(TBD.slice(0));
    TkDT.copy_size(TDT.slice(0));
    
    // Let's have option for parallelization
    int myEStarts = 0;
    int myEEnds = NE-1;
    for(int iE = myEStarts; iE <= myEEnds; ++iE){
        // Initialize the T(E) vector for energy E.
        uword Nk = mkMesh(iE).n_rows;
        mTek(iE).set_size(Nk);
        mTek(iE).zeros();
        double EE = mE(iE);

        // k-points assignment to CPUs
        int myKStarts, myKends;
        assignCPUs(Nk, myKStarts, myKends);
        
        for(int ik = myKStarts; ik <= myKends; ++ik){
            //---- calculate H(k) ----
            rowvec ki = mkMesh(iE).row(ik);
            // for the device
            //<botLayer|H|topLayer>
            generateHk(HkD, HD, ki, rd);
            // for bot contact
            generateHk(HkB, HB, ki, rd);
            // for top contact
            generateHk(HkT, HT, ki, rd);
            // <Bottom|H|Device>
            generateHk(TkBD, TBD, ki, rd);
            //<Device|H|Top>
            generateHk(TkDT, TDT, ki, rd);
            
            // compute transmission, T(k,E)
            mTek(iE)(ik) = computeTekNonRGF(HkB, TkBD, HkD, TkDT, HkT, EE);            
        }
#ifdef USE_MPI   
        // The master collects data ans progress report
        // and shows the porgress bar. 
        if(mIamMaster){
            // Collect T(E)
            mpi::reduce(*mWorld, mTek(iE), mTek(iE), plus<mat>(), mMasterId);
            
            // save T(E,k) as we progress
            if (iE%2 == 0){
                saveTek();
            }
            
            // Progress report
            ++progress;
        }else{
            // Send T(E) to the master
            mpi::reduce(*mWorld, mTek(iE), plus<mat>(), mMasterId);
        }
#else
        ++progress;
#endif

    }
    
    // DONE
    if(mIamMaster){
        progress.complete();
    }    
}

double Transport::computeTekNonRGF(const cx_mat& HkB, const cx_mat& TkBD, 
        const cx_mat& HkD, const cx_mat& TkDT, const cx_mat& HkT, double E)
{
    // ---- top contact surf green function ----
    uword nrPerLyrTpLead = HkT.n_rows/2;
    uvec lyr0TpLead = linspace<uvec>(0,nrPerLyrTpLead-1,nrPerLyrTpLead); 
    uvec lyr1TpLead = lyr0TpLead + nrPerLyrTpLead;
    cx_mat Tk01Tp = zeros<cx_mat>(HkT.n_rows, HkT.n_cols); // <supCell0|H|supCell1> for top lead
    Tk01Tp(lyr1TpLead,lyr0TpLead) = HkT(lyr1TpLead, lyr0TpLead); // <lyr1|H|lyr0> for top lead
    
    cx_mat SkT = eye<cx_mat>(HkT.n_rows, HkT.n_cols);
    cx_mat gst(HkT.n_rows, HkT.n_cols);
    if (!computeSurfG(gst, E, HkT, SkT, Tk01Tp, miEta, SurfGTolX)){
        cerr << "ERROR: Decimation did not converge, E = " << E << endl;
    }
    SkT.clear();
    Tk01Tp.clear();
    
    // ---- bottom contact surf green function ----
    uword nrPerLyrBtLead = HkB.n_rows/2;
    uvec lyr0BtLead = linspace<uvec>(0,nrPerLyrBtLead-1,nrPerLyrBtLead); 
    uvec lyr1BtLead = lyr0BtLead + nrPerLyrBtLead;
    cx_mat Tk10Bt = zeros<cx_mat>(HkB.n_rows, HkB.n_cols); // <supCell1|H|supCel0> for bottom lead
    Tk10Bt(lyr0BtLead,lyr1BtLead) = HkB(lyr0BtLead, lyr1BtLead); // <lyr0|H|lyr1> for bottom lead
    
    cx_mat SkB = eye<cx_mat>(HkB.n_rows, HkB.n_cols);
    cx_mat gsb(HkB.n_rows, HkB.n_cols);
    if (!computeSurfG(gsb, E, HkB, SkB, Tk10Bt, miEta, SurfGTolX)){
        cerr << "ERROR: Decimation did not converge, E = " << E << endl;
    }
    SkB.clear();
    Tk10Bt.clear();
    
    // ---- G11 ---- 
    cx_mat sig11Tp = TkDT*gst*trans(TkDT);
    gst.clear();
    
    // gamma11 bottom layer
    cx_mat sig11Bt = trans(TkBD)*gsb*TkBD;
    gsb.clear();
        
    // G11
    cx_mat G11 = E*eye<cx_mat>(HkD.n_rows, HkD.n_cols); // G11 = E*I
    G11 = inv(G11 - HkD - sig11Tp - sig11Bt); // G11 = inv[EI - H - SigL - SigR] 
    
    // T(E,k)
    sig11Bt = i*(sig11Bt - trans(sig11Bt)); // Gamma11
    sig11Tp = i*(sig11Tp - trans(sig11Tp)); // Gamma11
    cx_mat TEk = sig11Bt*G11*sig11Tp*trans(G11);
    
    return real(trace(TEk));
       
}


bool Transport::saveTek(){
    
    // create directory if not exist
    int status = mkdir(mOutputPath.c_str(), S_IRWXU | S_IXGRP);
    
    if(status != 0 && errno != EEXIST){
        throw runtime_error("Cannot create/access directory " + mOutputPath);
    }
    
    // save to a file
    string fileName = mOutputPath + "/" + mOutputFileName 
                      + "_" + Bias() + "." + mOutExt;
    ofstream transFile(fileName.c_str());
    if (!transFile.is_open()){
        throw runtime_error(" Failed to open file " + fileName + ".");
    }
    
    for (int iE = 0; iE < mTek.n_elem; ++ iE){
        if (!mTek(iE).empty()){
            transFile << "E = " << mE(iE) << endl;
            mat tek = mkMesh(iE);
            tek.insert_cols(2,mTek(iE));
            transFile << tek << endl;
        }
    }
    
    transFile.close();
    
    return true;
}

void Transport::readOptions(const TwOpts& opts){
    Calculator::readOptions(opts);
    TwOpts::CalculationOpts copts = opts.CalculationOptions();
    
    // Device type
    mDeviceType = copts.DeviceType;
            
    // get eta
    miEta.real() = 0;
    miEta.imag() = stod(copts.Eta);
    
    // inter layer distance
    mdocc = stod(opts.TightBindingParams().do0cc);
        
    double Emin = stod(copts.E.Start);
    double Emax = stod(copts.E.End);
    double delE = stod(copts.E.DelE);
    
    uword nE = 1;
    if (abs(delE) > 0){
        nE = uword(abs(Emax-Emin)/delE) + 1;
    }    
    mE = linspace<rowvec>(Emin, Emax, nE);
    
    // get the potential energy
    generatePotential(copts.Bias, copts.potential);
        
    // get the k-points
    createKMesh(copts.kpts);
    
    // Layer info
    extractLayerInfo(copts);
    
}

void Transport::createKMesh(const TwOpts::KptsOpts& kpts){
        
    int it;
    mat k;
    
    for (it = 0; it < kpts.KP.size(); ++it){
        genKPoints(kpts.KP[it], k);
    }
    
    for (it = 0; it < kpts.kl.size(); ++it){
        genKLines(kpts.kl[it], k);
    }
    
    for (it = 0; it < kpts.ks.size(); ++it){
        genKSurface(kpts.ks[it], k);
    }

    for (it = 0; it < kpts.kpar.size(); ++it){
        genKParallelo(kpts.kpar[it], k);
    }

    for (it = 0; it < kpts.ksc.size(); ++it){
        genKSmartCircle(kpts.ksc[it], k);
    }    
    
    for (it = 0; it < kpts.kc.size(); ++it){
        genKCircle(kpts.kc[it], k);
    }    
    
    mnKpoints = 0;
    // If KPoints, KLines, KSurfaces or KCircles are specified
    // then all energy points uses the same k-mesh. 
    if (!k.is_empty()){
        mkMesh.set_size(mE.n_cols);
        for (it = 0; it < mkMesh.n_elem; ++it){
            mkMesh(it) = k;
            mnKpoints += k.n_rows;
        }
    // If k-points
    // are auto generated from E-k then we use different k-mesh
    // for each energy.
    }else if (!kpts.ekc.isEmpty() || !kpts.eksc.isEmpty()){        
        
        string fileName;
        double centerX, centerY, searchSpan, delRad;
        if (!kpts.ekc.isEmpty()){
            fileName = kpts.ekc.FileName;
            // Dirac cone center
            centerX = stod(kpts.ekc.CenterX);
            centerY = stod(kpts.ekc.CenterY);
            searchSpan = stod(kpts.ekc.SearchSpan);
            delRad = stod(kpts.ekc.DelRadius);
            
        }else if(!kpts.eksc.isEmpty()){
            fileName = kpts.eksc.FileName;
            // Dirac cone center
            centerX = stod(kpts.eksc.CenterX);
            centerY = stod(kpts.eksc.CenterY);
            searchSpan = stod(kpts.eksc.SearchSpan);
            delRad = stod(kpts.eksc.DelK);
        }
        
        // Load the E-K
        mat ek;
        // only the master reads the file
        if(mIamMaster){
            ek.load(fileName + "_" + Bias() + "." + mOutExt);
        }   
#ifdef USE_MPI
        mpi::broadcast(*mWorld, ek, mMasterId);
#endif

        // extract the k-points and eigen energies
        rowvec kr = sqrt(pow(ek.row(X)-centerX, 2) + pow(ek.row(Y)-centerY,2));
        mat eigenE = ek.rows(2,ek.n_rows-1);
        ek.reset();
        
        double delE = searchSpan*miEta.imag();
        double delK = (300.0/1973.0)*delE;
        if(delK < delRad){
            delK = 5*delRad;
        }
        mkMesh.set_size(mE.n_cols);
        
        for (int iE = 0; iE < mE.n_cols; ++iE){
            double E = mE[iE];   
            mat newk;
            
            // Gather all the k-points that have energy near E
            uvec energyCircleIndx;        
            for(int iB = 0; iB < eigenE.n_rows; ++iB){
                rowvec bandE = eigenE.row(iB);
                uvec ec = find(abs(bandE-E) < delE);
                if(ec.empty()){
                    continue;
                }
                energyCircleIndx.insert_rows(energyCircleIndx.n_rows, ec);
            }
            if(energyCircleIndx.empty()){
                continue;
            }

            // Create a histogram of kr to see how many iso surfaces are 
            // crossing E
            vec circleK = kr.elem(energyCircleIndx);
            energyCircleIndx.clear();
            
            double minK = min(circleK);
            minK = minK < 0.0 ? 0.0 : minK;
            double maxK = max(circleK);
            
            uword nBeans = uword((maxK - minK)/delK) + 1;
            vec beans = linspace<vec>(minK, maxK, nBeans);
            uvec kfreq = hist(circleK, beans);
            
            // Lets find out how many k-grid zone we need          
            // Add a zero so that the following loop ends
            // correctly.
            bool inSlab = false; 
            vector<string> radiusStart;
            vector<string> radiusEnd;
            kfreq.insert_rows(kfreq.n_rows,1,true);
            for (int ib = 0; ib < kfreq.n_rows; ++ib){
                if (kfreq[ib] > 0){
                    if(!inSlab){
                        radiusStart.push_back(dtos(beans[ib] - delK));
                        inSlab = true;
                    }
                }else{
                    if(inSlab){
                        radiusEnd.push_back(dtos(beans[ib-1] + delK));
                        inSlab = false;
                    }
                }
            }

            if (radiusStart.size() != radiusEnd.size()){
                throw runtime_error("Cannot generate automatic k-grid.");
            }
            
            // Now create smart k-grid for transport calculation
            for (int is = 0; is < radiusStart.size(); ++is){
                if (!kpts.ekc.isEmpty()){
                    TwOpts::KCircle  kc;
                    kc.CenterX = kpts.ekc.CenterX;
                    kc.CenterY = kpts.ekc.CenterY;
                    kc.DelRadius = kpts.ekc.DelRadius;
                    kc.ThetaStart = kpts.ekc.ThetaStart;
                    kc.ThetaEnd = kpts.ekc.ThetaEnd;    
                    kc.DelTheta = kpts.ekc.DelTheta;
                    
                    kc.RadiusStart = radiusStart[is];
                    kc.RadiusEnd = radiusEnd[is];
                    
                    genKCircle(kc, newk);
                }else if (!kpts.eksc.isEmpty()){
                    TwOpts::KSmartCircle  ksc;
                    ksc.CenterX = kpts.eksc.CenterX;
                    ksc.CenterY = kpts.eksc.CenterY;
                    ksc.DelK = kpts.eksc.DelK;
                    ksc.ThetaStart = kpts.eksc.ThetaStart;
                    ksc.ThetaEnd = kpts.eksc.ThetaEnd;    

                    ksc.RadiusStart = radiusStart[is];
                    ksc.RadiusEnd = radiusEnd[is];
                    
                    genKSmartCircle(ksc, newk);
                }
            }


            mkMesh(iE) = newk;
            mnKpoints += newk.n_rows;
            radiusStart.clear();
            radiusEnd.clear();
            newk.reset();         
        }

    // If no k-points are generated then abort.    
    }else{
        throw runtime_error("No valid k-point specification found.");
    }
}

void Transport::extractLayerInfo(const TwOpts::CalculationOpts& copts){
    // Layer info
    stringstream ssbot(trim(copts.as.BotContact));
    stringstream ssdev(trim(copts.as.Device));
    stringstream sstop(trim(copts.as.TopContact));
    uword topStart, topEnd, botStart, botEnd;
    uword devStart, devEnd;
    
    // convert string to number
    ssbot >> botStart >> botEnd;
    sstop >> topStart >> topEnd;
    ssdev >> devStart >> devEnd;
    
    if (ssbot.fail() || sstop.fail() || ssdev.fail()){
        throw runtime_error("Contact and device information not fount.");
    }
    
    mBottomContact = linspace<uvec>(botStart, botEnd, botEnd - botStart + 1) - 1;
    mDevice = linspace<uvec>(devStart, devEnd, devEnd - devStart + 1) - 1;
    mTopContact = linspace<uvec>(topStart, topEnd, topEnd - topStart + 1) - 1;   
}

void Transport::genHam(cube& HB, cube& TBD, cube& HD, cube& TDT, cube& HT,
    vector<svec>& r){
    
    
    // For debugging lets save the contacts
    string neighFileName = mNeighFileName;
    
    // separate top and bottom layers
    Atoms supCell = *mSuperCell;
    Atoms topContact = supCell(mTopContact);
    Atoms botContact = supCell(mBottomContact);
    Atoms device = supCell(mDevice);
    
    mNeighFileName = "topLead_" + neighFileName;
    generateNeighHam(HT, r, &topContact);
    
    mNeighFileName = "botLead_" + neighFileName;
    generateNeighHam(HB, r, &botContact);
    
    mNeighFileName = "botLeadToDev_" + neighFileName;
    generateNeighHam(TBD, r, &botContact, &device);

    mNeighFileName = "devToTopLead_" + neighFileName;
    generateNeighHam(TDT, r, &device, &topContact);
    
    
    mNeighFileName = neighFileName;
    generateNeighHam(HD, r, &device);    
    
}




