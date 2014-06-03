/* 
 * File:   BandStruct.h
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on April 6, 2013, 9:25 AM
 * 
 * Description: Calculates band structure.
 * 
 */

#ifndef BANDSTRUCT_H
#define	BANDSTRUCT_H

#include "parallel/parloop.h"

#include "utils/ConsoleProgressBar.h"
#include "utils/Printable.hpp"
#include "utils/myenums.hpp"
#include "utils/NullDeleter.hpp"

#include "grid/grid.hpp"
#include "string/stringutils.h"
#include "maths/constants.h"
#include "maths/svec.h"
#include "maths/arma.hpp"

#include "atoms/Lattice.h"

#include <boost/smart_ptr.hpp>


#include <sys/types.h>
#include <sys/stat.h>



namespace qmicad{
namespace band{

using namespace maths::armadillo;
using namespace maths::constants;
using namespace maths::spvec;
using namespace utils;
using namespace utils::enums;
using namespace atoms;

namespace mpi = boost::mpi;
using boost::shared_ptr;
using std::string;

/**
 * Band structure parameters
 */
struct BandStructParams: public Printable{
    
    field<shared_ptr<cxmat> >H; //!< Hamiltonian of all nearest neighboring cells.
                                //!< H(0) is the cell # 0.
    field<shared_ptr<cxmat> >S; //!< Overlap matrix of all nearest neighboring cells.
                                //!< H(0) is the cell # 0.
    field<lcoord>   lc;         //!< Position vectors of neighbors.
    uint            nb;         //!< Number of bands to be saved.
    uint            no;         //!< Number of orbitals.
    uint            ne;         //!< Number of electrons in the super cell.
    uint            nn;         //!< Number of neighbors
    lvec            lv;         //!< Lattice vector.
    bool            isOrthogonal; //!< Is this orthogonal basis set?
    
    BandStructParams(uint nn, const string &prefix = ""):
        Printable(" " + prefix), H(nn), S(nn), lc(nn), nn(nn)
    {
        mTitle = "Band structure parameters";
    }
    
    void setH(shared_ptr<cxmat> H, uint it){ this->H(it) = H; }
    void setS(shared_ptr<cxmat> S, uint it){ this->S(it) = S; }
    void setLattCoord(shared_ptr<lcoord> pv, uint it) { this->lc(it) = *pv; }
};


/**
 * Band structure calculator.
 */
class BandStruct: public ParLoop {
// Fields    
protected:
    BandStructParams    mp;         //!< Simulation parameters.
    shared_ptr<mat>     mk;         //!< k-points:     (# of kpoints)  x  3.
    
    mat                 mE;         //!< Eigen energy: (# of kpoints)  x  (# of bands).
    mat                 mThisE;     //!< Eigen energy for this process: 
                                    //!< (# of kpoints for this proc)  x  (# of bands).
    int                 mlb;        //!< lowest band to calculate
    int                 mub;        //!< highest band to calculate    
    bool                mSaveAscii; //!< Save ASCII/Binary file?
    bool                mCalcEigV;  //!< Calculate eigen values?


    ConsoleProgressBar  mbar;//!< Progress bar.
// Methods    
public:
    BandStruct(shared_ptr<mat> pk, const BandStructParams &bp, 
            const Workers &workers, bool saveAscii = true);
   
    int     NumOfKpoints() const { return mN; };
    void    enableEigVec() { mCalcEigV = true; };
    void    save(string fileName);    

protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();

private:

};

}
}
#endif	/* BANDSTRUCTCALC_H */

