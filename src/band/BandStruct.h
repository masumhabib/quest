/* 
 * File:   BandStruct.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 6, 2013, 9:25 AM
 * 
 * Description: Calculates band structure.
 * 
 */

#ifndef BANDSTRUCT_H
#define	BANDSTRUCT_H

#include "../parallel/parloop.h"

#include "../utils/ConsoleProgressBar.h"
#include "../utils/Printable.hpp"
#include "../utils/myenums.hpp"
#include "../grid/grid.hpp"
#include "../string/stringutils.h"
#include "../maths/constants.h"
#include "../maths/svec.h"
#include "../maths/arma.hpp"

#include <boost/smart_ptr.hpp>


#include <sys/types.h>
#include <sys/stat.h>



namespace qmicad{
using namespace maths::armadillo;
using namespace maths::constants;
using namespace maths::spvec;
using namespace utils;
using namespace myenums;

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
    field<svec>     pv;         //!< Position vectors of neighbors.
    uint            Nb;         //!< Number of bands to be saved.
    uint            Ne;         //!< Number of electrons in the super cell.
    bool            saveAscii;  //!< Save ASCII/Binary file?
    
    void setH(shared_ptr<cxmat> H, uint it){ this->H(it) = H; }
    void setS(shared_ptr<cxmat> S, uint it){ this->S(it) = S; }
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


    ConsoleProgressBar  mbar;//!< Progress bar.
// Methods    
public:
    BandStruct(shared_ptr<mat> pk, const BandStructParams &bp, const Workers &workers);

    
    int NumOfKpoints() const { return mN; };
    
    
    // functionality
    bool save();    

protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    //virtual void    gather(vector<negfresult> &thisR, NegfResultList &all);

private:

};

}
#endif	/* BANDSTRUCTCALC_H */

