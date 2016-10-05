/* 
 * File:   computegs.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 22, 2014, 10:42 AM
 */

#ifndef COMPUTEGS_H
#define	COMPUTEGS_H

#include "maths/arma.hpp"
#include "maths/constants.h"

namespace quest{
namespace negf{

using namespace maths::armadillo;
using namespace maths::constants;

/** 
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


bool computegs(cxmat& gs, double E, const cxmat& Hii, const cxmat& Sii, 
        const cxmat& Tij, dcmplx ieta, double TolX);
}
}
#endif	/* COMPUTEGS_H */

