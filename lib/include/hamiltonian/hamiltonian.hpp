/* 
 * File:   hamiltonian.hpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 6, 2013, 5:52 PM
 * 
 * Description: Tight binding calculation logic and data.
 * 
 */

#ifndef HAMILTONIAN_HPP
#define	HAMILTONIAN_HPP

#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

#include "maths/constants.h"
#include "maths/arma.hpp"
#include "atoms/AtomicStruct.h"
#include "utils/vout.h"
#include "utils/std.hpp"

namespace quest{
namespace hamiltonian{

using namespace utils::stds;
using utils::Printable;
using namespace maths::armadillo;
using atoms::AtomicStruct;
using atoms::PeriodicTable;
using namespace maths::spvec;
using namespace maths::constants;
using utils::stds::static_pointer_cast;


template<class T>
class HamParams: public Printable{    
public:    
    HamParams(const string &prefix = ""):
        Printable(" " + prefix)
    {        
    }
    
    void   dtol(double dtol){ mdtol = dtol; update(); }
    double dtol() const { return mdtol; } 
    
    void   orthogonal(bool orth){ mortho = orth; update(); }
    bool   orthogonal() const { return mortho; }
    
    void   periodicTable(const PeriodicTable &pt) { mpt = pt; update(); }
    const  PeriodicTable &periodicTable() const { return mpt; }
    
    void   Bz(double Bz, int gauge = coord::X){ mBz = Bz; mBzGauge = gauge; update(); };
    double Bz() const { return mBz; }
    
    //!< Generate Hamiltonian between two atoms.
    virtual T twoAtomHam(const AtomicStruct& atomi, 
                            const AtomicStruct& atomj) const { return T(); };
    //!< Generate Overlap matrix between two atoms.
    virtual T twoAtomOvl(const AtomicStruct& atomi, 
                            const AtomicStruct& atomj) const { return T(); };
    
protected:
    // Updates internal parameters. Call it after changing any of the 
    // public parameters.
    virtual void update() {}; 

    // Caluclates Peierl's phase factor for magnetic field
    dcmplx calcPeierlsPhase(double xi, double yi, double xj, double yj) const;
    
protected:
    // Parameters required for all Hamiltonin
    double mdtol;         //!< Distance tolerance. 
    bool   mortho;        //!< Is this an orthogonal basis.
    PeriodicTable mpt;    //!< The periodic table required for this system.
    double mBz;           //!< The z-component of magnetic field.
    int    mBzGauge;      //!< gauge choice for the z-component.
    //!< pre-factor for magnetic phase = i*hbar/q/2
    double mfactor = 1E-20*q/hbar/2; 
    static constexpr double mBzTol = 1E-10;
};


// Caluclates Peierl's phase factor for magnetic field
template<class T>
dcmplx HamParams<T>::calcPeierlsPhase(double xi, double yi, 
        double xj, double yj) const 
{
    dcmplx phase = dcmplx(1,0);
    if(abs(mBz) > mBzTol) {
        dcmplx iphi = dcmplx(0,0);
        if (mBzGauge == coord::X) { // for A = (-Bz*y, 0, 0)
            iphi = dcmplx(0, mfactor*mBz*(xi - xj)*(yi + yj));
        } else if (mBzGauge == coord::Y) { // for A = (0, Bz*x, 0)
            iphi = dcmplx(0, mfactor*mBz*(yj - yi)*(xi + xj));
        }
        phase = exp(iphi);
    }
    return phase;
}

//!< Generates the hamiltonaina and overlap matrices between atomc block
//!< i and atomic block j. 
template<class T>
void generateHamOvl(T &hmat, T&smat, const HamParams<T> &p, 
        const AtomicStruct &bi, const AtomicStruct &bj){
    // Just for easy reference
    int nai = bi.NumOfAtoms();
    int naj = bj.NumOfAtoms();
    int noi = bi.NumOfOrbitals();
    int noj = bj.NumOfOrbitals();

    // Most of the matrix elements are zeros. So, we'll only change the 
    // non zero elements below.
    hmat = zeros<T>(noi, noj);
    smat = zeros<T>(noi, noj);

    // Lets find the neighbors. 
    int io = 0;
    for(int ia = 0; ia != nai; ++ia){
        AtomicStruct atomi = bi(ia);        // extract atom ia
        int ni = atomi.NumOfOrbitals();     // number of orbitals in atom ia

        int jo = 0;
        for(int ja = 0; ja != naj; ++ja){
            AtomicStruct atomj = bj(ja);        // extract atom ja    
            int nj = atomj.NumOfOrbitals();     // number of orbitals in atom ja
            // generate Hamiltonian matrix between orbitals of
            // atom i and atom j
            hmat(span(io,io+ni-1), span(jo,jo+nj-1)) = p.twoAtomHam(atomi, atomj);
            smat(span(io,io+ni-1), span(jo,jo+nj-1)) = p.twoAtomOvl(atomi, atomj);
            jo += nj;
        }

        io += ni;
    }    
}   

typedef HamParams<cxmat>  cxhamparams;
typedef HamParams<mat>    hamparams;

}
}
#endif	/* HAMILTONIAN_HPP */

