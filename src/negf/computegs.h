/* 
 * File:   computegs.h
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on January 22, 2014, 10:42 AM
 */

#ifndef COMPUTEGS_H
#define	COMPUTEGS_H

#include <armadillo>
#include "../utils/mymath.h"


using namespace arma;
using namespace constants;

bool computegs(cx_mat& gs, double E, const cx_mat& Hii, const cx_mat& Sii, 
        const cx_mat& Tij, dcmplx ieta, double TolX);

#endif	/* COMPUTEGS_H */

