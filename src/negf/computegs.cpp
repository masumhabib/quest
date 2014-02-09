/* 
 * File:   computegs.cpp
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on January 22, 2014, 10:42 AM
 */

#include "computegs.h"

namespace qmicad{
    
/**
 * 
 */
bool computegs(cxmat& gs, double E, const cxmat& Hii, const cxmat& Sii, 
        const cxmat& Tij, dcmplx ieta, double TolX){
 
    // initial guess (see the line just after Eq. 11 of [1])
    cxmat epi_1 = Hii;
    cxmat epsi_1 = Hii;
    cxmat alpai_1 = Tij;
    cxmat betai_1 = trans(Tij);
    double con_error_alpa = 10;
    double con_error_beta = 10;
    int iter = 0;
    
    bool  flag = true;

    // working matrices
    cxmat inv_mat, alpai, betai, epi, epsi;

    while((con_error_alpa > TolX) || (con_error_beta > TolX)){
        // ---- Eq. B6 of [1] (alpa == A & beta == B)---------
        inv_mat = inv((E+ieta)*Sii-epi_1);
        alpai = alpai_1*inv_mat*alpai_1;
        betai = betai_1*inv_mat*betai_1;
        epi = epi_1+alpai_1*inv_mat*betai_1+betai_1*inv_mat*alpai_1;
        epsi = epsi_1+alpai_1*inv_mat*betai_1;
        // ---- convergence checking (line 3rd after Eq. B6 of [1]) ---
        con_error_alpa = max(max(abs(alpai)));
        con_error_beta = max(max(abs(betai)));
        // ---- cycling variables ------
        alpai_1 = alpai;
        betai_1 = betai;
        epi_1 = epi;
        epsi_1 = epsi;
        // --- successful or not -----
        if (++iter > 500){
            flag = false;
            break;
        }
    }
    // ---- surface Green functions (Eq. B7 of [1])
    gs = inv(E*Sii-epsi);

    return flag;
}

}
