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

#include "parallel/Workers.h"

#include "utils/ConsoleProgressBar.h"
#include "utils/Printable.hpp"
#include "utils/myenums.hpp"
#include "utils/NullDeleter.hpp"
#include "utils/std.hpp"

#include "maths/grid.hpp"
#include "maths/constants.h"
#include "maths/svec.h"
#include "maths/arma.hpp"

#include "atoms/Lattice.h"

#include <sys/types.h>
#include <sys/stat.h>



namespace qmicad{
namespace band{

using namespace maths::armadillo;
using namespace maths::constants;
using namespace maths::spvec;
using namespace utils;
using namespace qmicad::parallel;
using namespace utils::enums;
using namespace atoms;

namespace mpi = boost::mpi;
using std::shared_ptr;
using std::string;

///**
// * Band structure parameters
// */
//struct BandStructParams: public Printable{
//    
//    
//    BandStructParams(uint nn, const string &prefix = ""):
//        Printable(" " + prefix), H(nn), S(nn), lc(nn), nn(nn)
//    {
//        mTitle = "Band structure parameters";
//    }
//    
//    void setH(shared_ptr<cxmat> H, uint it){ this->H(it) = H; }
//    void setS(shared_ptr<cxmat> S, uint it){ this->S(it) = S; }
//    void setLattCoord(shared_ptr<lcoord> pv, uint it) { this->lc(it) = *pv; }
//};


/**
 * Band structure calculator.
 */
class BandStruct: public Printable {
// Methods    
public:
    BandStruct(const Workers &workers, uint nn, bool orthoBasis = true, 
            bool calcEigV = false, const string &prefix = "");
   
    void    nb(uint nb);
    uint    nb(){ return mnb; };
    void    ne(uint ne);
    uint    ne(){ return mne; };
    void    lv(const lvec& lv);
    lvec    lv() { return mlv; }
    void    lc(field<lcoord> lc);
    void    lc(const lcoord &lc, int ineigh);
    void    k(const mat &k);
    mat     k(){ return mk; };
    
    void    H(const field<shared_ptr<cxmat> > &H);
    void    S(const field<shared_ptr<cxmat> > &S);    
    void    H(shared_ptr<cxmat> H, int ineigh);
    void    S(shared_ptr<cxmat> S, int ineigh);
    
    int     NumOfKpoints() const { return mN; };
    void    enableEigVec() { mCalcEigV = true; };
    
    void    run();
    void    save(string fileName, bool saveAsText = true);    

protected:
    void    resetBandIndices();
    
    virtual void    prepare();
    virtual void    preCompute(long il);
    virtual void    compute(long il);  
    virtual void    postCompute(long il);
    virtual void    collect();
    

private:

// Fields    
protected:
    const Workers       &mWorkers;  //!< MPI worker processes.
    long                mN;         //!< Number of grid points.
    long                mMyN;       //!< Number of grid points to be calculated by this process.    
    long                mMyStart;   //!< Start point for this CPU.
    long                mMyEnd;     //!< End point for this CPU.

    uint                mnn;        //!< Number of neighbors
    uint                mnb;        //!< Number of bands to be saved.
    uint                mno;        //!< Number of orbitals.
    uint                mne;        //!< Number of electrons in the super cell.
    bool                mOrthoBasis;//!< Is this orthogonal basis set?
    field<lcoord>       mlc;        //!< Position vectors of neighbors.
    lvec                mlv;        //!< Lattice vector.

    mat                 mk;         //!< k-points:     (# of kpoints)  x  3.    
    field<shared_ptr<cxmat> >mH;    //!< Hamiltonian of all nearest neighboring cells.
                                    //!< H(0) is the cell # 0.
    field<shared_ptr<cxmat> >mS;    //!< Overlap matrix of all nearest neighboring cells.
                                    //!< H(0) is the cell # 0.    
    
    mat                 mE;         //!< Eigen energy: (# of kpoints)  x  (# of bands).
    mat                 mThisE;     //!< Eigen energy for this process: 
                                    //!< (# of kpoints for this proc)  x  (# of bands).
    long                mlb;        //!< lowest band to calculate
    long                mub;        //!< highest band to calculate    
    bool                mSaveAscii; //!< Save ASCII/Binary file?
    bool                mCalcEigV;  //!< Calculate eigen values?


    ConsoleProgressBar  mbar;//!< Progress bar.

};

}
}
#endif	/* BANDSTRUCTCALC_H */

