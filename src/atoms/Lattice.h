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

#include "../maths/svec.h"
#include "../utils/std.hpp"
#include "../utils/serialize.hpp"
#include "../utils/Printable.hpp"

namespace qmicad{
using namespace utils::stds;
using utils::Printable;
using namespace maths::armadillo;
using namespace maths::spvec;

/* Lattice coordinates */
struct LatticeCoordinate: public Printable{
    //fields
    int n1;
    int n2;
    int n3;

    LatticeCoordinate(int n1 = 0, int n2 = 0, int n3 = 0, const string &prefix = "");

    //operators
    LatticeCoordinate& operator+= (const LatticeCoordinate& rhs);
    friend LatticeCoordinate operator+ (LatticeCoordinate lhs, const LatticeCoordinate& rhs);
    LatticeCoordinate& operator-= (const LatticeCoordinate& rhs);
    friend LatticeCoordinate operator- (LatticeCoordinate lhs, const LatticeCoordinate& rhs);
    
    virtual string toString() const; 
};

/* Lattice vector */
struct LatticeVector: public Printable {
    // three vectors. format: a = a0 x + a1 y + a2 z
    svec a1;    // Vector 1
    svec a2;    // Vector 2
    svec a3;    // Vector 3
    
    LatticeVector(const string &prefix = "");
    void zeros();
    
    double la1(){return sqrt(sum(square(a1))); }
    double la2(){return sqrt(sum(square(a2))); }
    double la3(){return sqrt(sum(square(a3))); }
    
    //operators
    LatticeVector& operator+= (const LatticeVector& rhs);
    friend LatticeVector operator+ (LatticeVector lhs, const LatticeVector& rhs);
    LatticeVector& operator-= (const LatticeVector& rhs);
    friend LatticeVector operator- (LatticeVector lhs, const LatticeVector& rhs);    
    
    friend svec operator* (const LatticeVector &lv, const LatticeCoordinate& lc);    
    
    virtual string toString() const; 
    
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

