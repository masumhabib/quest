/* 
 * File:   CohRgfa.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 * 
 * Created on April 22, 2013, 8:05 PM
 */

#include "CohRgfa.h"

namespace qmicad{
CohRgfa::CohRgfa(NegfParams newp, double E, string newprefix):
        Printable(newprefix), mp(newp), mE(E), 
        mN(newp.nb-2), miLc(0), miRc(newp.nb-1),
        mDi(this, miLc, miRc, newp.DCache), 
        mTl(this, miLc, miRc+1, newp.TCache),
        mgrc(this, miLc+1, miRc, newp.grcCache),
        mglc(this, miLc, miRc-1, newp.glcCache),
        mGii(this, miLc+1, miRc-1, newp.GiiCache),
        mGi1(this, miLc+2, miRc-1, newp.Gi1Cache),
        mGiN(this, miLc+1, miRc-2, newp.GiNCache),
        mGiip1(this, miLc+1, miRc-2, newp.Giip1Cache),
        mGiim1(this, miLc+2, miRc-1, newp.Giim1Cache)
{
    mTitle = "Transport";
    // Fermi functions
    mf0 = fermi(mE, mp.muS, mp.kT);
    mfNp1 = fermi(mE, mp.muD, mp.kT);
}

CohRgfa::~CohRgfa() {
}

/*
 * Transmission operator T(E) = tr{Gamma_1,1*[A_1,1 - G_1,1*Gamma_1,1*G_1,1']}
 */
cxmat CohRgfa::TEop(uint traceOverN){
    // Full Green function G_1,1
    // G_1,1 = [D_1,1 - sig_l_1,1 - T_1,2*grc_2,2*T_2,1]^-1
    // sig_l_1,1 = T_1,0*glc_0,0*T_0,1
    const cxmat &glc00 = mglc(0);                    // Get or caluclate glc_0,0
    cxmat Gaml11 = mTl(1)*glc00*trans(mTl(1)); 
    Gaml11 = i*(Gaml11 - trans(Gaml11));
    const cxmat &G11 = mGii(1);                       // Get or caluclate G_1,1
    
    cxmat TEop = Gaml11*(i*(G11 - trans(G11)) - G11*Gaml11*trans(G11));    
    return trace<cxmat>(TEop, traceOverN);
}

/*
 * Electron density operator. It returns Gn_i,i or sum_j(Gn_i,j*Sj,i)
 * depending on the orthogonality of the basis set.
 */
cxmat CohRgfa::niOp(uint traceOverN){
    
}

/*
 * Current operator at terminal 1.
 * I1op
 */
cxmat CohRgfa::I1Op(uint traceOverN){
    
    // Full Green function G_1,1
    // G_1,1 = [D_1,1 - sig_l_1,1 - T_1,2*grc_2,2*T_2,1]^-1
    // sig_l_1,1 = T_1,0*glc_0,0*T_0,1
    const cxmat &glc00 = mglc(0);                    // Get or caluclate glc_0,0
    cxmat Sigl11 = mTl(1)*glc00*trans(mTl(1)); 
    const cxmat &G11 = mGii(1);                       // Get or caluclate G_1,1
    
    cxmat Gn11;
    cxmat I1op;
    // Density matrix: Gn11 = G^n_1,1
    // G^n_1,1 = Al_1,1*(f1-fN) + [A_1,1]*fN
    // Al_1,1 = G_1,1*gamma_1,1*G1,1'
    // A_1,1 = i*(G_1,1 - G_1,1')
    // gamma_1,1 =  i*(Sigl_1,1-Sigl_1,1')    
    cxmat Gaml11 = i*(Sigl11 - trans(Sigl11)); 
    Gn11 = G11*Gaml11*trans(G11)*(mf0-mfNp1) + i*(G11 - trans(G11))*mfNp1;
    // Current operator
    I1op = Gn11*trans(Sigl11) - Sigl11*Gn11 
          +G11*Gaml11*mf0 - Gaml11*trans(G11)*mf0;

    return trace<cxmat>(i*I1op, traceOverN);    
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
        NegfParams &p = mnegf->mp;
        double E = mnegf->mE;
        double VL = (*p.V(iLc))(0); // all the atoms on a contact have the save bias
        computegs(glci, E+VL, *p.H0(iLc), *p.S0(iLc), Tiim1, p.ieta, p.SurfGTolX);
    // calculate glc_i,i using recursive equation:
    // glc_i = [ES_ii - H_ii - U_ii - T_ii-1*glc_i-1*T_i-1i]^-1;
    }else{
        glci = inv(mnegf->mDi(ib) - Tiim1*glcim1*trans(Tiim1));
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
        NegfParams &p = mnegf->mp;
        double E = mnegf->mE;
        double VR = (*p.V(iRc))(0); // all the atoms on a contact have the save bias
        computegs(grci, E+VR, *p.H0(iRc), *p.S0(iRc), trans(Tip1i), p.ieta, 
                  p.SurfGTolX);

    // Calculate grc_i,i using recursive equation:
    // grc_i = [ES_ii - H_ii - U_ii - T_ii+1*grc_i+1*T_i+1i]^-1;
    }else{
        grci = inv(mnegf->mDi(ib) - trans(Tip1i)*grcip1*Tip1i);        
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
    cxmat& Hl = *(mnegf->mp.Hl(ii));
    
    // for orthogonal basis
    if (mnegf->mp.isOrthogonal){
        Tl = Hl;
    // for non-orthogonal basisi
    }else{
        cxmat& Sl = *(mnegf->mp.Sl(ii));
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
    cxmat &Hii = *(mnegf->mp.H0(ii));
    
    // for orthogonal basis
    if (mnegf->mp.isOrthogonal){
        cxmat EIii = eye<cxmat>(Hii.n_rows, Hii.n_cols);
        vec &Vii = *(mnegf->mp.V(ii));
        EIii.diag().fill(E);    // E*I
        cxmat UIii(Hii.n_rows, Hii.n_cols, fill::zeros);
        UIii.diag() = dcmplx(-1, 0)*Vii;
        Dii = EIii - UIii - Hii;
        
    // for non-orthogonal basis
    }else{
        cxmat& Sii = *(mnegf->mp.S0(ii));
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
    cxmat &Sl = *(mp.Sl(i));
    Ul.copy_size( Sl);
    for(int m = 0; m < Ul.n_rows; ++m){
        for(int n = 0; n < Ul.n_cols; ++n){
            // if we are at the contacts: U_0,-1 and U_N+2,N+1
            // then use potential of the contacts V(0) and V(N+1) respectively.
            if (i == miLc || i == miRc+1){
                vec &Vi = *(mp.V(i));
                Ul(m, n) = -(Vi(m) + Vi(n))/2* Sl(m,n);
            }else{
                vec &Vi = *(mp.V(i));
                vec &Vim1 = *(mp.V(i-1));
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
    cxmat &Si = *(mp.S0(i));
    Ui.copy_size(Si);
    for(int m = 0; m < Ui.n_rows; ++m){
        for(int n = 0; n < Ui.n_cols; ++n){
            //[Uij]m,n = - (V_im+V_in)/2*[Sij]_m,n
            vec &Vi = *(mp.V(i));
            Ui(m, n) = -(Vi(m) + Vi(n))/2*Si(m,n);
        }
    }
    return Ui;
}

}


