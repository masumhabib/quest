/* 
 * File:   Transport.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 22, 2013, 8:05 PM
 */

#ifndef TRANSPORT_H
#define	TRANSPORT_H

#include <sys/types.h>
#include <sys/stat.h>

#include <string>
#include <armadillo>
#include <vector>

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

#include "Calculator.h"
#include "../tightBinding/TightBinder.h"
#include "../atoms/Atoms.h"
#include "../utils/ConsoleProgressBar.h"
#include "../utils/mymath.h"
#include "../utils/stringutils.h"
#include "../utils/computeSurfG.h"

using namespace std;
using namespace arma;
using namespace constants;

#ifdef USE_MPI
namespace mpi = boost::mpi;
#endif


class Transport: public Calculator {
// Fields    
protected:
    dcmplx miEta;
    
    rowvec mE;          // 1       x (# of energy points)
    // The field mkMesh stores k-points for each energy point 
    // in a (# of k-points) x 2 matrix where the first row stores kx and the 
    // last row stores ky
    field<mat> mkMesh;  // (# of energy points) matrices of (# of k-points) x 2 size
    // The field mTek stores k-points and transmission for each energy
    // in a (# of k-points) x 1 vector.
    // row stores transmission.
    field<colvec> mTek; // (# of energy points) matrices of (# of k-points) x 1 size
    
    // Layer indices
    ucolvec mTopContact;
    ucolvec mDevice;
    ucolvec mBottomContact;
    
    uword mnKpoints;
    
    string mDeviceType;
    
    double mdocc;
    
    static const double SurfGTolX = 1E-8;
        
public:
#ifdef USE_MPI
    Transport(mpi::communicator *world, Atoms* atm, TightBinder *tb, string prefix = "");
#else
    Transport(Atoms* atm, TightBinder *tb, string prefix = "");
#endif
    virtual ~Transport();
    
    //Access functions
    uword NumOfEnergyPoints() const { return mE.n_cols; };
    uword NumOfKPoints() const { return mnKpoints; };
    
    // functionality
    void computeTek();
    void computeTekNonRGF();
    bool saveTek();
    void readOptions(const TwOpts& opts);
protected:
        // Helpers
        double computeTekNonRGF(const cx_mat& HkB, const cx_mat& TkBD, 
               const cx_mat& HkD, const cx_mat& TkDT, const cx_mat& HkT, 
               double E);
        void createKMesh(const TwOpts::KptsOpts& kpts);
        void extractLayerInfo(const TwOpts::CalculationOpts& copts);
        
        void genHam(cube& HB, cube& TBD, cube& HD, cube& TDT, cube& HT,
                    vector<svec>& r);
private:
    Transport(); // We need atoms an tight binder from the start

};

#endif	/* TRANSPORT_H */

