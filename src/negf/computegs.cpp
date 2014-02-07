/* 
 * File:   computegs.cpp
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on January 22, 2014, 10:42 AM
 */

#include "computegs.h"

/* 
 * This function calculates the surface green's function of a
 * semi-infinite structure using decimation technique.
 *=====================================================================
 *
 * [1] M. Galperin, et al., "Numerical Calculation of Tunneling Fluxes," 
 *     J. Chem. Phys., vol. 117, pp. 10817-10826, 2002.
 *
 *======================================================================
 * E ---------> energy in eV
 * Hii -------> Hamiltonian of the principal layer
 * Sii -------> overlap matrix of the orbitals of principal layer
 * Tij -------> coupling matrix between 0th and 1st principal layers of
 *              right semi-infinite lead (Tij = tij-ESij)
 * ieta ------> imaginary energy (eV)
 * TolX ------> convergence tolerance factor
 * gs --------> surface green function (by default gs = gR; for gs = gL 
 *              send Tij dagar)
 * returns ---> false if not converged in 500 iteration
 *====================================================================
*/

bool computegs(cx_mat& gs, double E, const cx_mat& Hii, const cx_mat& Sii, 
        const cx_mat& Tij, dcmplx ieta, double TolX){
 
    // initial guess (see the line just after Eq. 11 of [1])
    cx_mat epi_1 = Hii;
    cx_mat epsi_1 = Hii;
    cx_mat alpai_1 = Tij;
    cx_mat betai_1 = trans(Tij);
    double con_error_alpa = 10;
    double con_error_beta = 10;
    int iter = 0;
    
    bool  flag = true;

    // working matrices
    cx_mat inv_mat, alpai, betai, epi, epsi;

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
