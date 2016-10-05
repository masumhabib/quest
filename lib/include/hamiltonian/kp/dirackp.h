/* 
 * File:   dirackp.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 10:16 AM
 */

#ifndef DIRACKP_H
#define	DIRACKP_H

#include <armadillo>

#include "maths/constants.h"
#include "atoms/AtomicStruct.h"
#include "utils/Printable.hpp"
#include "hamiltonian/hamiltonian.hpp"

namespace quest{
namespace hamiltonian{

using namespace maths::armadillo;
using namespace maths::constants;

/**
 * Dirac k.p Hamiltonian parameters.
 * See: Phys. Rev. B 86, 085131 (2012) and the references therein.
 */
class DiracKpParams: public cxhamparams {
public:
    DiracKpParams(const string &prefix = "");
    
    double a() const {return ma; }
    void   a(double newa ){ ma = newa; update(); }
    double K() const {return mK; }
    void   K(double newK ){ mK = newK; update(); }

    //!< Generate Hamiltonian between two atoms.
    virtual cxmat twoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj) const;    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat twoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj) const;
protected:
    //!< Updates internal tight binding parameters calculated using 
    //!< k.p model. Call it after changing any of the k.p parameters.
    virtual void update() {};
    
private:
 

protected:
    //!< k.p parameters
    double ma;            //!< discretization length in x direction
    double mK;            //!< parameter to avoid fermion doubling

    //!< calculated tight binding matrices
    cxmat mI;
    cxmat meps;
    cxmat mt01x;
    cxmat mt10x;
    cxmat mt01y;
    cxmat mt10y;
    
    
};

}
}
#endif	/* DIRACKP_H */

