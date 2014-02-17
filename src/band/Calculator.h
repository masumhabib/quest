/* 
 * File:   Calculator.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 23, 2013, 11:21 AM
 * 
 * Description: Base of all calculator classes.
 * 
 */

#ifndef CALCULATOR_H
#define	CALCULATOR_H

#include <sys/types.h>
#include <sys/stat.h>

#include <armadillo>
#include <vector>
#include <string>
#include <stdexcept>

#include <boost/mpi.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

#include "../atoms/Atoms.h"
#include "../options/TwOpts.h"
#include "../tightBinding/TightBinder.h"
#include "../utils/mymath.h"
#include "../utils/Printable.hpp"

using namespace std;
using namespace arma;
using namespace constants;

#ifdef USE_MPI
namespace mpi = boost::mpi;
#endif

class Calculator: public Printable {
// Fields
protected:
    //options 
    static const string mOutExt;
    string  mOutputFileName;
    string  mOutputPath;
    string  mNeighFileName;
    bool    mIsSaveNeigh;      // generate and save neighboring supercells?

#ifdef USE_MPI
    static const int    PROGRESS = 0;
    static const int    mMasterId = 0;
    mpi::communicator*  mWorld;
#endif
    int                 mnCpu;
    int                 mMyCpuId;
    bool                mIamMaster;
    
    Atoms*       mSuperCell;
    TightBinder* mTightBinder;  
    double       mV;             // bias voltage
    colvec       mPotential;  // potential


// Methods
public:
#ifdef USE_MPI
    Calculator(mpi::communicator *world, Atoms* atm, TightBinder *tb, string prefix = "");
#else
    Calculator(Atoms* atm, TightBinder *tb, string prefix = "");
#endif
    virtual ~Calculator();
    
    //Access
    string Bias() const;
    
    // operations
    virtual string toString() const;
    void readOptions(const TwOpts& opts);
    void assignCPUs(int N, int& myStartPoint, int& myEndPoint);
    
protected:
    void generateNeighHam(cube &H, vector<svec>& r, 
            const Atoms* superCellI = NULL, const Atoms* superCellJ = NULL);
    void generateHk(cx_mat& Hk, cube& H, rowvec& k, vector<svec>& r);
    void generatePotential(const string& bias, const vector<TwOpts::Potential>& potential);

    void genKPoints(const string& kp, mat& k);
    void genKLines(const TwOpts::KLine& kl, mat& k);
    void genKSurface(const TwOpts::KSurface& ks, mat& k);
    void genKParallelo(const TwOpts::KParallelo& kpar, mat& k);
    void genKSmartCircle(const TwOpts::KSmartCircle& ksc, mat& k);
    void genKCircle(const TwOpts::KCircle& kc, mat& k);
    
    
private:
    Calculator(); // dont want a default constructor

};

#endif	/* CALCULATOR_H */

