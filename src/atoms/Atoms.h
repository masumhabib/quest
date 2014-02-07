/* 
 * File:   Atoms.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 5, 2013, 1:27 PM
 * 
 * Description: The structure Atom represents an atom in the periodic
 * table. The class Atoms encapsulates everything about the atoms present in
 * the device.
 * 
 */

#ifndef ATOMS_H
#define	ATOMS_H

#include <iostream>
#include <fstream>
#include <vector>
#include <iomanip> 
#include <armadillo>

#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include "Lattice.h"
#include "../utils/stringutils.h"
#include "../utils/svec.h"


using namespace std;
using namespace arma;
using spacevec::svec;
    
    
        
// Atom in the periodic table
struct Atom {
    int ia;     // atomic number
    string sym; // atomic symbol
    int ne;     // number of electrons
    int no;     // number of orbitals

    Atom() {};
    Atom(int ia, string sym, int ne, int no){
        this->ia = ia;
        this->sym = sym;
        this->ne = ne;
        this->no = no;
    };
    Atom(const Atom& orig):
    ia(orig.ia),
    sym(orig.sym),
    ne(orig.ne),
    no(orig.no){
        
    };
    
// For MPI send/receive
private:
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & ia;
        ar & sym;
        ar & ne;
        ar & no;
    }

};


/* All atoms in the structure */
class Atoms {
// Fields
protected:
    vector<Atom> mPeriodicTable;  // our periodic table

    int mNumAtoms;                // no of atoms
    int mNumOrbitals;             // no of orbitals
    int mNumElectrons;            // no of electrons

    // atoms and their coordinates
    icolvec mAtomId;            // id in periodic table
    mat     mXyz;        // atom coordinates
    lvec    mLatticeVector;     // lattice vector

// Methods    
public:
    // Constructors and destructors.
    Atoms(); // default
    // constructs from a GaussView gjf file
    Atoms(const string& gjfFileName ); 
    // constructs from a GaussView gjf file and a periodic table
    Atoms(const string& gjfFileName, const vector<Atom>& PeriodicTable);
    // constructs from a atomic coordinates and lattice vector
    Atoms(const icolvec& atomId, const mat& coordinate, const lvec& lv);
    // constructs from a atomic coordinates and lattice vector
    Atoms(const icolvec& atomId, const mat& coordinate, const lvec& lv,
    const vector<Atom>& PeriodicTable);
    // copy constructor
    Atoms(const Atoms& orig);
    virtual ~Atoms(){};
    friend void swap(Atoms& first, Atoms& second);
    // operators
    Atoms  operator()(span s) const;                  // Get a sub cell 
    Atoms  operator()(const ucolvec& index) const;    // Get a sub cell
    Atoms  operator()(uword i) const;                 // Get one atom cell
    
    Atoms& operator= (Atoms rhs);                     // assignment
    Atoms& operator+= (const Atoms& atoms);           // concatenation
    Atoms& operator-= (const lcoord& latticeCoord);   // coordinate shifting
    Atoms& operator-= (const svec& positionVect);     // coordinate shifting
    Atoms& operator+= (const lcoord& latticeCoord);   // coordinate shifting
    Atoms& operator+= (const svec& positionVect);     // coordinate shifting

    // non-mamber operators
    friend Atoms operator- (Atoms atm, const lcoord& latticeCoord);     // atm2 = atm1 - lc;
    friend Atoms operator- (Atoms atm, const svec& positionVect);       // atm2 = atm1 - r;
    friend Atoms operator+ (Atoms atm, const lcoord& latticeCoord);     // atm2 = atm1 + lc;
    friend Atoms operator+ (Atoms atm, const svec& positionVect);       // atm2 = atm1 + r;
    friend Atoms operator+ (Atoms atmi, const Atoms& atmj);             // concatanation
    friend ostream& operator<< (ostream& out, const Atoms &b);

    // utilities
    void importGjf(const string &gjfFileName);
    void exportGjf(const string &gjfFileName);

    // access functions
    // get's
    Atom        AtomAt(uword i) const;       // Get one atom at i
    string      Symbol(int i) const { return convertIndexToSym(mAtomId(i)); } ;
    double      X(int i) const { return mXyz(i, spacevec::X); };
    double      Y(int i) const { return mXyz(i, spacevec::Y); };
    double      Z(int i) const { return mXyz(i, spacevec::Z); };
    vec         X() const { return mXyz.col(spacevec::X); };
    vec         Y() const { return mXyz.col(spacevec::Y); };
    vec         Z() const { return mXyz.col(spacevec::Z); };
    int         NumOfAtoms() const { return mNumAtoms; };
    int         NumOfOrbitals() const { return mNumOrbitals; };
    int         NumOfElectrons() const { return mNumElectrons; };
    const lvec& LatticeVector() const { return mLatticeVector; };
    //set's
    void        LatticeVector(const lvec& a) { this->mLatticeVector = a; };
    void        PeriodicTable(const vector<Atom> &periodicTable);


protected:
    void        init();
    void        initPeriodicTable();
    int         computeNumOfOrbitals();
    int         computeNumOfElectrons();

private:
    int         convertSymToIndex(const string& sym) const;
    string      convertIndexToSym(int ind) const;

    // For MPI send/receive
    friend class boost::serialization::access;

    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & mPeriodicTable;
        ar & mNumAtoms;
        ar & mNumOrbitals;
        ar & mNumElectrons;
        ar & mAtomId;
        ar & mXyz;
        ar & mLatticeVector;
    }
    
};

#endif	/* ATOMS_H */

