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

#include <boost/smart_ptr.hpp>

#include <string>
#include <armadillo>
#include <vector>

#include <sys/types.h>
#include <sys/stat.h>



namespace qmicad{
using namespace maths::armadillo;
using namespace maths::constants;
using namespace maths::spacevec;
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
    field<svec>     pv;         //!< Position vectors of neighbors.
    uint            Nb;         //!< Number of bands to be saved.
    bool            saveAscii;  //!< Save ASCII/Binary file?
};


/**
 * Band structure calculator.
 */
class BandStruct: public ParLoop {
// Fields    
protected:
    BandStructParams    mp;  //!< Simulation parameters.
    shared_ptr<mat>     mk;  //!< k-points: (# of kpoints) x        2 (or 3).
    mat                 mE;  //!< Eigen energy: (# of bands)  x  (# of kpoints).

    ConsoleProgressBar  mbar;//!< Progress bar.
// Methods    
public:
    BandStruct(shared_ptr<mat> pk, const BandStructParams &bp, const Workers &workers);

    
    int NumOfKpoints() const { return mk->n_rows; };
    
    
    // functionality
    bool save();    

protected:
    virtual void    prepare();
    virtual void    preCompute(int il);
    virtual void    compute(int il);  
    virtual void    postCompute(int il);
    virtual void    collect();
    virtual void    stepCompleted();
    //virtual void    gather(vector<negfresult> &thisR, NegfResultList &all);

private:

};

}
#endif	/* BANDSTRUCTCALC_H */

