/* 
 * File:   CohRgfa.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 22, 2013, 8:05 PM
 */

#include "negf/CohRgfa.h"

namespace quest{
namespace negf{

CohRgfa::CohRgfa(uint nb, double kT, dcmplx ieta, bool orthogonal, string newprefix):
        Printable(newprefix),
        mnb(nb), mkT(kT), mieta(ieta), morthogonal(orthogonal),
        mH0(nb), mS0(nb), mHl(nb+1), mSl(nb+1), mV(nb),
        mN(nb-2), miLc(0), miRc(nb-1),
        mDi(this, miLc, miRc), 
        mTl(this, miLc, miRc+1),
        mgrc(this, miLc+1, miRc),
        mglc(this, miLc, miRc-1),
        mGii(this, miLc+1, miRc-1),
        mGi1(this, miLc+2, miRc-1),
        mGiN(this, miLc+1, miRc-2),
        mGiip1(this, miLc+1, miRc-2),
        mGiim1(this, miLc+2, miRc-1)
{
    mTitle = "Coherent Transport using RGF";
}

// set chemical potential
void CohRgfa::mu(double muD, double muS){
    mmuS = muS;
    mmuD = muD;
    mf0 = fermi(mE, mmuS, mkT);
    mfNp1 = fermi(mE, mmuD, mkT);
}


void CohRgfa::E(double E){
    mE = E;
    reset();

    mf0 = fermi(mE, mmuS, mkT);
    mfNp1 = fermi(mE, mmuD, mkT);    
}


void CohRgfa::H(const field<shared_ptr<cxmat> > &H0, const field<shared_ptr<cxmat> > &Hl){
    if (H0.n_elem != mnb){
        throw runtime_error("In CohRgfa::H(): size of H0 should be equal to number of blocks.");
    }

    if (Hl.n_elem != mnb+1){
        throw runtime_error("In CohRgfa::H(): size of Hl should be equal to number of blocks + 1.");
    }
    
    mH0 = H0;
    mHl = Hl;    
}

void CohRgfa::S(const field<shared_ptr<cxmat> > &S0, const field<shared_ptr<cxmat> > &Sl){
    if (S0.n_elem != mnb){
        throw runtime_error("In CohRgfa::S(): size of S0 should be equal to number of blocks.");
    }

    if (Sl.n_elem != mnb+1){
        throw runtime_error("In CohRgfa::S(): size of Sl should be equal to number of blocks + 1.");
    }

    mS0 = S0;
    mSl = Sl;
}

void CohRgfa::V(const field<shared_ptr<vec> > &V){
    if (V.n_elem != mnb){
        throw runtime_error("In CohRgfa::V(): size of V should be equal to number of blocks.");
    }
    
    mV = V;
}

string CohRgfa::toString() const {
    stringstream out;
    out << Printable::toString() << ":" << endl;
    out << mPrefix << " IsOrthogonal = " << (morthogonal ? "Yes" : "No")  << endl;
    out << mPrefix << " nb           = " << mnb << endl;
    out << mPrefix << " N            = " << mN << endl;
    out << mPrefix << " ieta         = " << mieta << endl;
    out << mPrefix << " kT           = " << mkT << endl;
    out << mPrefix << " muS          = " << mmuS << endl;
    out << mPrefix << " muD          = " << mmuD;

    return out.str();
}

/*
 * Hole density operator. It returns A_ii - Gn_i,i or sum_j(A_i,j*Sj,i - Gn_i,j*Sj,i)
 * depending on the orthogonality of the basis set.
 */
cxmat CohRgfa::pOp(uint N, int ib, ucol *atomsTracedOver){
    cxmat pOp(N, N, fill::zeros);
    // If out of range, sum over the device
    if (ib >= miRc || ib <= miLc){
        for (uint ib = miLc+1; ib < miRc; ++ib){
            pOp += trace<cxmat>(Aop(ib, ib) - Gn(ib, ib), N, atomsTracedOver);
        }
    }else{
        pOp = trace<cxmat>(Aop(ib, ib) - Gn(ib, ib), N, atomsTracedOver);
    }
    
    return pOp/(2*pi);    
}

/*
 * Electron density operator. It returns 
 * Gn_i,i or sum_j(Gn_i,j*Sj,i) depending on the orthogonality of the basis set.
 */
cxmat CohRgfa::nOp(uint N, int ib, ucol *traveOveratoms){

    cxmat nOp(N, N, fill::zeros);
    // If out of range, sum over the device
    if (ib >= miRc || ib <= miLc){
        for (uint ib = miLc+1; ib < miRc; ++ib){
            nOp += trace<cxmat>(Gn(ib, ib), N, traveOveratoms);
        }
    }else{
        nOp = trace<cxmat>(Gn(ib, ib), N, traveOveratoms);
    }
    
    return nOp/(2*pi);        
}

/*
 * Density of States. 
 */
cxmat CohRgfa::DOSop(uint N, ucol *atomsTracedOver){
    cxmat D(N, N, fill::zeros);
    for (uint ib = miLc+1; ib < miRc; ++ib){
        D += Aop(N, ib, atomsTracedOver);
    }  
    return D/(2*pi);
}

/*
 * Spectral function for block ib
 */
cxmat CohRgfa::Aop(uint N, uint ib, ucol *traveOveratoms){
    cxmat A = mGii(ib);
    A = i*(A - trans(A));
    
    return trace<cxmat>(A, N, traveOveratoms);    
}


/*
 * Generic current operators: current from block i to block j.
 * -----------------------------------------------------------------------------
 */
cxmat CohRgfa::Iop(uint N, uint ib, uint jb, ucol *atomsTracedOver){
        
    if (ib == miLc && jb == miLc + 1){
        // current from contact 0 to device.
        return i*trace<cxmat>(I0op(), N, atomsTracedOver);
//        return I0op(N, atoms);
    }else if (ib == miLc + 1 && jb == miLc){
        // current from device to contact 0.
        return -i*trace<cxmat>(I0op(), N, atomsTracedOver);
//        return -I0op(N, atoms);
    } else if (ib == miRc && jb == miRc-1){
        // current from contact N to device.
        return i*trace<cxmat>(INop(), N, atomsTracedOver);
//        return INop(N, atoms);
    } else if (ib == miRc-1 && jb == miRc){
        // current from device to contact N.
        return -i*trace<cxmat>(INop(), N, atomsTracedOver);
//        return -INop(N, atoms);
    } else if (abs(int(ib) - int (jb)) <= 1 && ib > miLc && ib < miRc && jb > miLc && jb < miRc) {
        return i*trace<cxmat>(Iijop(ib, jb), N, atomsTracedOver);
//        return Iijop(N, ib, jb, atoms);
    }else{
        stringstream error;
        error << "ERROR: Iop(i, j) only works for |i-j| <= 1 and  0<= i,j <= N+1."
              << endl << "But we've got, i = " << ib << ", j = " << jb 
              << " and N = " << N << ".";
        throw invalid_argument(error.str());        
    }
}

/*
 * Transmission operator T(E) = tr{Gamma_1,1*[A_1,1 - G_1,1*Gamma_1,1*G_1,1']}
 * -----------------------------------------------------------------------------
 */
cxmat CohRgfa::TEop(uint N, ucol *traveOveratoms){
    // Full Green function G_1,1
    // G_1,1 = [D_1,1 - sig_l_1,1 - T_1,2*grc_2,2*T_2,1]^-1
    // G_1,1 = [D_1,1 - sig_l_1,1 - SigL_1,1]^-1    
    const cxmat &G11 = mGii(1);                       // Get or caluclate G_1,1
    const cxmat &Gaml11 = GamL11();
    cxmat G11a = trans(G11);
    cxmat TEop = Gaml11*(i*(G11 - G11a) - G11*Gaml11*G11a); 
    return trace<cxmat>(TEop, N, traveOveratoms);
}

/*
 * Current from block # i and block j where i and j are neighbors.
 * -----------------------------------------------------------------------------
 */
inline cxmat CohRgfa::Iijop(uint ib, uint jb){
    cxmat Iijop;
    cxmat &Gnij = Iijop;
    
    Gnij = Gn(ib, jb);
    //I_i,j = H_i,j*Gn_j,i - Gn_i,j*H_j,i; 
    if (ib < jb){
        //I_i,i+1 = H_i,i+1*Gn_i+1,i - Gn_i,i+1*H_i+1,i
        Iijop = trans(mTl(jb))*trans(Gnij) - Gnij*mTl(jb);
    }else if (ib > jb){
        //I_i+1,i = H_i+1,i*Gn_i,i+1 - Gn_i+1,i*H_i,i+1
        Iijop = mTl(ib)*trans(Gnij) - Gnij*trans(mTl(ib));
    }else{
        Iijop = mDi(ib)*Gnij - Gnij*mDi(ib);
    }
    
    //return i*trace<cxmat>(Iijop, N);    
    return Iijop;
}


/*
 * Current operator at right contact (block # N).
 * INOp
 * -----------------------------------------------------------------------------
 */
inline cxmat CohRgfa::INop(){

    const cxmat &SigrNN = SigRNN();
    const cxmat &GamrNN = GamRNN();
    const cxmat &GNN = mGii(mN);            // Get or caluclate G_N,N
    cxmat GNNa = trans(GNN);
    
    // Density matrix: Gn11 = G^n_1,1
    // G^n_1,1 = Al_1,1*(f1-fN) + [A_1,1]*fN
    // Al_1,1 = G_1,1*gamma_1,1*G1,1'
    // A_1,1 = i*(G_1,1 - G_1,1')
    cxmat GnNN = GNN*GamrNN*GNNa*(mfNp1-mf0) + i*(GNN - GNNa)*mf0;
    // Current operator
    cxmat INop = GnNN*trans(SigrNN) - SigrNN*GnNN 
          + GNN*GamrNN*mfNp1 - GamrNN*GNNa*mfNp1;

    //return i*trace<cxmat>(INop, N, atoms);    
    return INop;
}


/*
 * Current operator at terminal 1.
 * I1op
 * -----------------------------------------------------------------------------
 */
inline cxmat CohRgfa::I0op(){
    
    const cxmat &Sigl11 = SigL11();
    const cxmat &Gaml11 = GamL11();
    const cxmat &G11 = mGii(1);            // Get or caluclate G_1,1
    cxmat G11a = trans(G11);
    
    // Density matrix: Gn11 = G^n_1,1
    // G^n_1,1 = Al_1,1*(f1-fN) + [A_1,1]*fN
    // Al_1,1 = G_1,1*gamma_1,1*G1,1'
    // A_1,1 = i*(G_1,1 - G_1,1')
    cxmat Gn11 = G11*Gaml11*G11a*(mf0-mfNp1) + i*(G11 - G11a)*mfNp1;
    // Current operator
    cxmat I0op = Gn11*trans(Sigl11) - Sigl11*Gn11 
          +G11*Gaml11*mf0 - Gaml11*G11a*mf0;
    
    //return i*trace<cxmat>(I1op, N, atoms);    
    return I0op;
}

/**
 * Correlation function.
 * ----------------------------------------------------------------------------- 
 */
inline cxmat CohRgfa::Gn(uint ib, uint jb){
    cxmat Gnij;

    // Gn_i,j = i*[G_i,j - G_j,i']*fN + G_i,1*Gam_1,1*G_j,1'*(f1-fN)    
    Gnij = (i*mfNp1)*(G(ib, jb) - trans(G(jb, ib))) 
              + (mf0 - mfNp1)*(G(ib,1)*GamL11()*trans(G(jb,1)));
       
    return Gnij;
}

/*
 * Full retarded green function.
 * -----------------------------------------------------------------------------
 */
inline const cxmat& CohRgfa::G(uint ib, uint jb){
    if (ib == jb){
        return mGii(ib);
    }else if (ib + 1 == jb){
        return mGiip1(ib);
    }else if (ib == jb + 1){
        return mGiim1(ib);
    }else if(jb == miLc + 1){
        return mGi1(ib);
    }else if (jb == miRc - 1){
        return mGiN(ib);
    }else{
        stringstream error;
        error << "ERROR: G( " << ib << ", " << jb << ") is not implemented";
        throw runtime_error(error.str());
    }
}

/*
 * The Giim1 class members to compute:
 * G_i,i-1 = grc_i,i*T_i,i-1*G_i-1,i-1
 * =============================================================================
 */

/* 
 * This function calculates full Green function G along the diagonal,
 * upper diagonal and lower diagonal from block iGii to block ib.
 * -----------------------------------------------------------------------------
 * ib --------> Block index for which we want G_i,i, G_i,i+1 and G_i,i-1.
 * -----------------------------------------------------------------------------
 */
const cxmat& CohRgfa::Giim1::operator ()(int ib){
    cxmat& Giim1 = getAt(ib);
    if (!isStored(ib)){
        computeGiim1(Giim1, mnegf->mGii(ib-1), ib);
    }
    return Giim1;
}

/* 
 * This function calculates the full Green function along the lower
 * diagonal: G_i,i-1
 * -----------------------------------------------------------------------------
 * Giim1 -----> Output: G_i,i-1
 * Gim1im1 ---> Input: G_i-1,i-1
 * ib --------> Block index for which we want G_i,i-1.
 * -----------------------------------------------------------------------------
 */
inline void CohRgfa::Giim1::computeGiim1(cxmat& Giim1, const cxmat& Gim1im1, int ib){
    CohRgfa &nf = *mnegf;
    // Calculate G_i,i-1 using recursive equation    
    // G_i,i-1 = grc_i,i*T_i,i-1*G_i-1,i-1
    Giim1 = nf.mgrc(ib)*nf.mTl(ib)*Gim1im1;
    mIt = ib;
}


/*
 * The Giip1 class members to compute:
 * G_i,i+1 = G_i,i*T_i,i+1*grc_i+1,i+1
 * =============================================================================
 */

/* 
 * This function calculates full Green function G along the upper diagonal
 * for block ib.
 * -----------------------------------------------------------------------------
 * ib --------> Block index for which we want G_i,i+1
 * -----------------------------------------------------------------------------
 */
const cxmat& CohRgfa::Giip1::operator ()(int ib){
    cxmat& Giip1 = getAt(ib);
    if (!isStored(ib)){
        computeGiip1(Giip1, mnegf->mGii(ib), ib);
    }
    return Giip1;
}

/* 
 * This function calculates the full Green function along the upper
 * diagonal: G_i,i+1
 * -----------------------------------------------------------------------------
 * Giip1 -----> Output: G_i,i+1
 * Gii -------> Input: G_i,i.
 * ib --------> Block index for which we want G_i,i+1.
 * -----------------------------------------------------------------------------
 */
inline void CohRgfa::Giip1::computeGiip1(cxmat& Giip1, const cxmat& Gii, int ib){
    CohRgfa &nf = *mnegf;
    // Calculate G_i,i+1 using recursive equation    
    // G_i,i+1 = G_i,i*T_i,i+1*grc_i+1,i+1
    Giip1 = Gii*trans(nf.mTl(ib+1))*nf.mgrc(ib+1);
    mIt = ib;
}


/*
 * The GiN class members to compute:
 * G_i,N = grc_i,i*T_i,i+1*G_i+1,N
 * G_N-1,N is calculated from G_N,N
 * =============================================================================
 */

/* 
 * This function calculates and returns full Green function G along
 * column#N from iGiN to ib.
 * -----------------------------------------------------------------------------
 * ib --------> Block index for which we want grc.
 * -----------------------------------------------------------------------------
 */
const cxmat& CohRgfa::GiN::operator ()(int ib){
    if (!isStored(ib)){
        // Block from which we start the calculation is the one just before
        // the last calculated block.
        int igStart = mIt - 1; 
        for (int ig = igStart; ig >= ib; --ig){
            int igp1 = (ig == end()) ? ig:ig+1;
            cxmat &GiN = getAt(ig);
            cxmat &Gip1N = getAt(igp1);
            computeGiN(GiN, Gip1N,ig);
        }
    }
    return getAt(ib);
}


/* 
 * This function calculates the full Green function along 
 * column N: G_i,N
 * -----------------------------------------------------------------------------
 * GiN -------> Output: G_i,N
 * Gip1N -----> Input: G_i+1,N
 * ib --------> Block index for which we want G_i,N.
 * -----------------------------------------------------------------------------
 */
inline void CohRgfa::GiN::computeGiN(cxmat& GiN, const cxmat& Gip1N, int ib){
    
    CohRgfa &nf = *mnegf;
    // Calculate GNm1N from GNN = Gii(N)
    if (ib == nf.mGii.end() - 1){ 
        GiN = nf.mgrc(ib)*trans(nf.mTl(ib+1))*nf.mGii(ib+1);
        
    // Calculate G_i,N using recursive equation    
    }else{
        // G_i,N = grc_i,i*T_i,i+1*G_i+1,N
        GiN = nf.mgrc(ib)*trans(nf.mTl(ib+1))*Gip1N;        
    }
    mIt = ib;    
}


/*
 * The Gi1 class members to compute:
 * G_i,1 = grc_i,i*T_i,i-1*G_i-1,1
 * G_2,1 is calculated from G_1,1
 * =============================================================================
 */

/* 
 * This function returns the full Green function G along
 * column#1 for block ib. If Gi1 was not previously calculated,
 * this function calculates it.
 * ib --------> Block index for which we want G_i,1.
 */
const cxmat& CohRgfa::Gi1::operator ()(int ib){
    if (!isStored(ib)){
        // Block from which we start the calculation is the one just after
        // the last calculated block.
        int igStart = mIt + 1; 
        //if (igStart < mIt){
        //    igStart = mnegf->miLc;
        //}
        for (int ig = igStart; ig <= ib; ++ig){
            int igm1 = (ig == begin()) ? ig:ig-1;
            cxmat &Gi1 = getAt(ig);
            cxmat &Gim11 = getAt(igm1);
            computeGi1(Gi1, Gim11,ig);
        }
    }
    return getAt(ib);
}

/* 
 * This function calculates the full Green function along 
 * column 1: G_i,1
 * Gi1 -------> Output: G_i,1
 * Gim11 -----> Input: G_i-1,1
 * ib --------> Block index for which we want G_i,1.
 */
inline void CohRgfa::Gi1::computeGi1(cxmat& Gi1, const cxmat& Gim11, int ib){
    CohRgfa &nf = *mnegf;
    // Calculate G21 from G11 = Gii(1)
    if (ib == nf.mGii.begin() + 1){ 
        Gi1 = nf.mgrc(ib)*nf.mTl(ib)*nf.mGii(ib-1);
        
    // Calculate G_i,1 using recursive equation        
    }else{
        // G_i,1 = grc_i,i*T_i,i-1*G_i-1,1
        Gi1 = nf.mgrc(ib)*nf.mTl(ib)*Gim11;
    }
    mIt = ib;    
}

/*
 * The Gii class members to compute:
 * G_1,1 = [ES_1,1 - H_1,1 - U_1,1 - sig1_1,1 - sig2_1,1]^-1
 * and
 * G_i,i = grc_i,i + grc_i,i*T_i,i-1*G_i-1,i-1*T_i-1,i*grc_i,i
 * =============================================================================
 */

/* 
 * This function returns full Green function G along the diagonal
 * for block # ib. If it was not previously calculated, this function
 * calculates it.
 * ib --------> Block index for which we want G_i,i.
 */
const cxmat& CohRgfa::Gii::operator ()(int ib){
    if (!isStored(ib)){
        // Block from which we start the calculation is the one just after
        // the last calculated block.
        int igStart = mIt + 1; 
        //if (igStart < mIt){
        //    igStart = mnegf->miLc;
        //}
        for (int ig = igStart; ig <= ib; ++ig){
            int igm1 = (ig == begin()) ? ig:ig-1;
            cxmat &Gii = getAt(ig);
            cxmat &Giim1 = getAt(igm1);
            computeGii(Gii, Giim1,ig);
        }
    }
    return getAt(ib);    
}

/* 
 * This function calculates the full Green function along the 
 * diagonal: G_i,i
 * Gii -------> Output: G_i,i
 * Gim1im1 ---> Input: G_i-1,i-1.
 * ib --------> Block index for which we want G_i,i.
 */
inline void CohRgfa::Gii::computeGii(cxmat& Gii, const cxmat& Gim1im1, int ib){
    // Calculated G_1,1 using
    // G_1,1 = [ES_1,1 - H_1,1 - U_1,1 - sig1_1,1 - sig2_1,1]^-1
    CohRgfa &nf = *mnegf;
    if(ib == nf.miLc+1){
        const cxmat &Tiim1 = nf.mTl(ib);
        const cxmat &Tip1i = nf.mTl(ib+1);
        Gii = inv(nf.mDi(ib) - Tiim1*nf.mglc(ib-1)*trans(Tiim1) 
                - trans(Tip1i)*nf.mgrc(ib+1)*Tip1i);
        
    // Otherwise,
    // calculate G_i,i using recursive equation    
    // G_i,i = grc_i,i + grc_i,i*T_i,i-1*G_i-1,i-1*T_i-1,i*grc_i,i
    }else{
        const cxmat &grci = nf.mgrc(ib);
        const cxmat &Tiim1 = nf.mTl(ib);
        Gii = grci + grci*Tiim1*Gim1im1*trans(Tiim1)*grci;
    }    
    mIt = ib;
}


/*
 * The glc class members to compute:
 * glc_i = [ES_ii - H_ii - U_ii - T_ii-1*glc_i-1*T_i-1i]^-1;
 * =============================================================================
 */
 
/* This function returns glc for block # it. If glc(it) was not stored 
 * in the memory from a previous calculation then it calculates 
 * all the left connected Green function from
 * last calculated block upto it.
 * it --------> Block index for which we want glc.
 */
const cxmat& CohRgfa::glc::operator ()(int ib){
    if (!isStored(ib)){
        // Block from which we start the calculation is the one just after
        // the last calculated block.
        int igStart = mIt + 1; 
        //if (igStart < mIt){
        //    igStart = mnegf->miLc;
        //}
        for (int ig = igStart; ig <= ib; ++ig){
            int igm1 = (ig == begin()) ? ig:ig-1;
            cxmat &glci = getAt(ig);
            cxmat &glcim1 = getAt(igm1);
            computeglc(glci, glcim1,ig);
        }
    }
    return getAt(ib);    
}

/* 
 * This function calculates the left connected Green function: glc
 * glci ------> Output: glc_i,i
 * glcim1 ------> Input: glc_i-1,i-1.
 * ib --------> Block index for which we want glc. If ib is the index of 
 *              any one of the contacts, this function will calculate
 *              surface Green function.
 */
inline void CohRgfa::glc::computeglc(cxmat& glci, const cxmat& glcim1, int ib){
    int iLc = mnegf->miLc;
    const cxmat &Tiim1 = mnegf->mTl(ib); //load T_ib,ib-1
    // If this is the left contact, calculate surface Green function.
    if(ib == iLc){
        double E = mnegf->mE;
        double VL = (*mnegf->mV(iLc))(0); // all the atoms on a contact have the save bias
        computegs(glci, E+VL, *mnegf->mH0(iLc), *mnegf->mS0(iLc), Tiim1, mnegf->mieta, CohRgfa::SurfGTolX);
    // calculate glc_i,i using recursive equation:
    // glc_i = [ES_ii - H_ii - U_ii - T_ii-1*glc_i-1*T_i-1i]^-1;
    // glc_i = [ES_ii - H_ii - U_ii - SigL_ii]^-1;
    }else{
        // Calculate or load sigma_1,1
        if (ib == (mnegf->miLc + 1)){
            // save sigma_1,1
            if (mnegf->mSigL11.empty()){
                mnegf->computeSigL(mnegf->mSigL11, Tiim1, glcim1);
            }
            glci = mnegf->mSigL11;
        // Calculate SigL_i,i
        }else{
            mnegf->computeSigL(glci, Tiim1, glcim1);
        }
        glci = inv(mnegf->mDi(ib) - glci);
    }    
    mIt = ib;
}

/* 
 * The grc class members.
 * grc_i = [ES_ii - H_ii - U_ii - T_ii+1*grc_i+1*T_i+1i]^-1;
 * =============================================================================
 */
 
/* This function returns grc for block # it. If grc(it) was not stored 
 * in the memory from a previous calculation then it calculates 
 * all the right connected Green function from
 * last calculated block upto it.
 * it --------> Block index for which we want grc.
 */
const cxmat& CohRgfa::grc::operator ()(int ib){
    if (!isStored(ib)){
        // Block from which we start the calculation is the one just before
        // the last calculated block.
        int igStart = mIt - 1; 
        //if (igStart > mIt){
        //    igStart = mnegf->miRc;
        //}
        for (int ig = igStart; ig >= ib; --ig){
            int igp1 = (ig == end()) ? ig:ig+1;
            cxmat &grci = getAt(ig);
            cxmat &grcip1 = getAt(igp1);
            computegrc(grci, grcip1,ig);
        }
    }
    return getAt(ib);    
}

/* 
 * This function calculates the right connected Green function: grc.
 * grci ------> Output: grc_i,i
 * grcip1 ----> Input: grc_i+1,i+1.
 * it --------> Block index for which we want grc. If it is the index of 
 *              any one of the contacts, this function will calculate
 *              surface Green function.
 */
inline void CohRgfa::grc::computegrc(cxmat& grci, const cxmat& grcip1, int ib){
    int iRc = mnegf->miRc; 
    const cxmat &Tip1i = mnegf->mTl(ib+1); //load T_ib+1,ib
    // If this is the right contact then calculate surface Green function.
    if(ib == iRc){
        double E = mnegf->mE;
        double VR = (*mnegf->mV(iRc))(0); // all the atoms on a contact have the save bias
        computegs(grci, E+VR, *mnegf->mH0(iRc), *mnegf->mS0(iRc), trans(Tip1i), mnegf->mieta, 
                  CohRgfa::SurfGTolX);

    // Calculate grc_i,i using recursive equation:
    // grc_i = [ES_ii - H_ii - U_ii - T_ii+1*grc_i+1*T_i+1i]^-1;
    // grc_i = [ES_ii - H_ii - U_ii - SigR_ii]^-1;
    }else{
        if (ib == mnegf->mN){
            // save sigma_N,N
            if (mnegf->mSigRNN.empty()){
                mnegf->computeSigR(mnegf->mSigRNN, Tip1i, grcip1);
            }
            grci = mnegf->mSigRNN;
        }else{
            mnegf->computeSigR(grci, Tip1i, grcip1);
        }
        grci = inv(mnegf->mDi(ib) - grci);        
    }
    mIt = ib;
}


/*
 * Tl class:
 * Tij = [Hij + USij - ESij]
 * =============================================================================
 */
const cxmat& CohRgfa::Tl::operator ()(int ib){
    cxmat& M = getAt(ib);
    if (!isStored(ib)){
        computeTl(M, ib);
    }
    return M;
}

inline void CohRgfa::Tl::computeTl(cxmat& Tl, int ib){
    int ii = toArrayIndx(ib);
    double E = mnegf->mE;
    cxmat& Hl = *(mnegf->mHl(ii));
    
    // for orthogonal basis
    if (mnegf->morthogonal){
        Tl = Hl;
    // for non-orthogonal basisi
    }else{
        cxmat& Sl = *(mnegf->mSl(ii));
        cxmat Ul = mnegf->Ul(ii);
        Tl = Hl + Ul - E*Sl;
    }
    mIt = ib;
}

/*
 * Di class:
 * Dii = [ESii - USii - Hii]
 * =============================================================================
 */
const cxmat& CohRgfa::Di::operator ()(int ib){
    
    cxmat& M = getAt(ib);
    if (!isStored(ib)){
        computeDi(M, ib);
    }
    return M;
}

inline void CohRgfa::Di::computeDi(cxmat& Dii, int ib){
    int ii = toArrayIndx(ib);
    double E = mnegf->mE;
    cxmat &Hii = *(mnegf->mH0(ii));
    
    // for orthogonal basis
    if (mnegf->morthogonal){
        cxmat EIii = eye<cxmat>(Hii.n_rows, Hii.n_cols);
        vec &Vii = *(mnegf->mV(ii));
        EIii.diag().fill(E);    // E*I
        cxmat UIii(Hii.n_rows, Hii.n_cols, fill::zeros);
        UIii.diag() = dcmplx(-1, 0)*Vii;
        Dii = EIii - UIii - Hii;
        
    // for non-orthogonal basis
    }else{
        cxmat& Sii = *(mnegf->mS0(ii));
        cxmat USii = mnegf->Ui(ii);
        Dii = E*Sii - USii - Hii;
    }
    mIt = ib;
}


/*
 * Lower diagonal of U matrix for non-orthogonal basis
 * [Uij]m,n = - (V_im+V_in)/2*[Sij]_m,n
 */
inline cxmat CohRgfa::Ul(int i){
    cxmat Ul;
    cxmat &Sl = *(mSl(i));
    Ul.copy_size( Sl);
    for(int m = 0; m < Ul.n_rows; ++m){
        for(int n = 0; n < Ul.n_cols; ++n){
            // if we are at the contacts: U_0,-1 and U_N+2,N+1
            // then use potential of the contacts V(0) and V(N+1) respectively.
            if (i == miLc || i == miRc+1){
                vec &Vi = *(mV(i));
                Ul(m, n) = -(Vi(m) + Vi(n))/2* Sl(m,n);
            }else{
                vec &Vi = *(mV(i));
                vec &Vim1 = *(mV(i-1));
                Ul(m, n) = -(Vi(m) + Vim1(n))/2* Sl(m,n);
            }
        }
    }
    return Ul;
}
/*
 * Diagonal blocks of U matrix for non-orthogonal basis
 */
inline cxmat CohRgfa::Ui(int i){
    cxmat Ui;
    cxmat &Si = *(mS0(i));
    Ui.copy_size(Si);
    for(int m = 0; m < Ui.n_rows; ++m){
        for(int n = 0; n < Ui.n_cols; ++n){
            //[Uij]m,n = - (V_im+V_in)/2*[Sij]_m,n
            vec &Vi = *(mV(i));
            Ui(m, n) = -(Vi(m) + Vi(n))/2*Si(m,n);
        }
    }
    return Ui;
}

/*
 * SigL_i,i = T_ii-1*glc_i-1*T_i-1i
 */
inline void CohRgfa::computeSigL(cxmat& SigLii, const cxmat& Tiim1, const cxmat& glcim1){
    SigLii = Tiim1*glcim1*trans(Tiim1);
}

/*
 * SigR_i,i = T_ii+1*grc_i+1*T_i+1i 
 */
inline void CohRgfa::computeSigR(cxmat& SigRii, const cxmat& Tip1i, const cxmat& grcip1){
    SigRii = trans(Tip1i)*grcip1*Tip1i;
}

/*
 * sigL_1,1 = T_1,0*glc_0,0*T_0,1
 */
inline const cxmat& CohRgfa::SigL11(){
    if (mSigL11.empty()){
        // sigL_1,1 = T_1,0*glc_0,0*T_0,1        
        computeSigL(mSigL11, mTl(1), mglc(0)); 
    }
    return mSigL11;
}

/*
 * SigR_N,N = T_N,N+1*grc_N+1,N+1*T_N+1,N 
 */
inline const cxmat& CohRgfa::SigRNN(){
    if (mSigRNN.empty()){
        computeSigR(mSigRNN, mTl(mN+1), mgrc(mN+1));
    }
    return mSigRNN;
}

inline const cxmat& CohRgfa::GamL11(){
    if (mGamL11.empty()){
        mGamL11 = i*(SigL11() - trans(SigL11()));
    }
    return mGamL11;
}

inline const cxmat& CohRgfa::GamRNN(){
    if (mGamRNN.empty()){
        mGamRNN = i*(SigRNN() - trans(SigRNN()));
    }
    return mGamRNN;
}

inline void CohRgfa::reset(){
    mDi.reset();
    mTl.reset();
    mgrc.reset();
    mglc.reset();
    mGii.reset();
    mGi1.reset();
    mGiN.reset();
    mGiip1.reset();
    mGiim1.reset();
    mSigL11.reset();
    mSigRNN.reset();
    mGamL11.reset();
    mGamRNN.reset();
}

}
}




