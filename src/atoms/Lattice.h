/* 
 * File:   Lattice.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 7, 2013, 9:46 AM
 * 
 * Description: Lattice vector definition
 * 
 */

#ifndef LATTICE_H
#define	LATTICE_H

#include <iostream>
#include <armadillo>
#include <boost/serialization/access.hpp>

#include "../maths/svec.h"
#include "../utils/serialize.hpp"

namespace qmicad{
using std::ostream;
using std::endl;
using namespace maths::armadillo;
using namespace maths::spvec;

/* Lattice coordinates */
struct LatticeCoordinate{
    //fields
    int n1;
    int n2;
    int n3;

    LatticeCoordinate();
    virtual ~LatticeCoordinate();
    LatticeCoordinate(int n1, int n2, int n3);

    //operators
    LatticeCoordinate& operator+= (const LatticeCoordinate& rhs);
    friend LatticeCoordinate operator+ (LatticeCoordinate lhs, const LatticeCoordinate& rhs);
    LatticeCoordinate& operator-= (const LatticeCoordinate& rhs);
    friend LatticeCoordinate operator- (LatticeCoordinate lhs, const LatticeCoordinate& rhs);
    
    friend ostream& operator << (ostream & out, const LatticeCoordinate & lc);
};

/* Lattice vector */
struct LatticeVector {
    // three vectors. format: a = a0 x + a1 y + a2 z
    svec a1;    // Vector 1
    svec a2;    // Vector 2
    svec a3;    // Vector 3
    
    LatticeVector();
    virtual ~LatticeVector();
    void zeros();
    
    //operators
    LatticeVector& operator+= (const LatticeVector& rhs);
    friend LatticeVector operator+ (LatticeVector lhs, const LatticeVector& rhs);
    LatticeVector& operator-= (const LatticeVector& rhs);
    friend LatticeVector operator- (LatticeVector lhs, const LatticeVector& rhs);    
    
    friend ostream& operator << (ostream & out, const LatticeVector& lv);
    
    // For MPI send/receive
    private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version)
    {
        ar & a1;
        ar & a2;
        ar & a3;
    }
};

typedef struct LatticeVector lvec;
typedef struct LatticeCoordinate lcoord;

}

#endif	/* LATTICE_H */

