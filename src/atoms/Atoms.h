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

#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include "../string/stringutils.h"
#include "../maths/svec.h"
#include "../maths/arma.hpp"
#include "../grid/grid.hpp"

#include "Lattice.h"

namespace qmicad{
using utils::trim;
using utils::MatGrid; 
using namespace utils::stds;
using namespace maths::armadillo;
namespace spacevec = maths::spacevec;
    
    
        
/**
 * The Atom structure describes properties of an atom as seen in 
 * the periodic table.
 */
struct Atom {
    int ia;     //!< Atomic number
    string sym; //!< Atomic symbol
    int ne;     //!< Number of electrons
    int no;     //!< Number of orbitals

    //! Default constructor.
    Atom() {};
    //! Detailed constructor taking all the properties.
    Atom(int ia, string sym, int ne, int no){
        this->ia = ia;
        this->sym = sym;
        this->ne = ne;
        this->no = no;
    };
    
    //! Copy constructor.
    Atom(const Atom& orig):
    ia(orig.ia),
    sym(orig.sym),
    ne(orig.ne),
    no(orig.no){
        
    };
    
private:
    //! For MPI send/receive.
    friend class boost::serialization::access;

    //! For MPI send/receive.
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & ia;
        ar & sym;
        ar & ne;
        ar & no;
    }

};


/** 
 * All atoms in the structure.
 */
class Atoms {
// Fields
protected:
    vector<Atom> mPeriodicTable;  // our periodic table

    int mNumAtoms;                // no of atoms
    int mNumOrbitals;             // no of orbitals
    int mNumElectrons;            // no of electrons

    // atoms and their coordinates
    icol mAtomId;            // id in periodic table
    mat     mXyz;        // atom coordinates
    lvec    mLatticeVector;     // lattice vector

// Methods    
public:
    //! Constructors and destructors.
    Atoms(); // default
    //! constructs from a GaussView GJF file.
    Atoms(const string& gjfFileName ); 
    //! constructs from a GaussView GJF file and a periodic table.
    Atoms(const string& gjfFileName, const vector<Atom>& PeriodicTable);
    //! Constructor for k.p fake atoms.
    Atoms(double Lx, double Ly, double ax, double ay, const vector<Atom>& periodicTable);
    //! Constructs from a atomic coordinates and lattice vector.
    Atoms(const icol& atomId, const mat& coordinate, const lvec& lv);
    //! Constructs from a atomic coordinates and lattice vector.
    Atoms(const icol& atomId, const mat& coordinate, const lvec& lv,
    const vector<Atom>& PeriodicTable);
    //! Copy constructor.
    Atoms(const Atoms& orig);
    virtual ~Atoms(){};
    friend void swap(Atoms& first, Atoms& second);
    // operators
    Atoms  operator()(span s) const;                  // Get a sub cell 
    Atoms  operator()(const ucol& index) const;    // Get a sub cell
    Atoms  operator()(uint i) const;                 // Get one atom cell
    
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
    //! Import atoms from Gaussview file.
    void importGjf(const string &gjfFileName);
    //! Export atoms to Gaussview file.
    void exportGjf(const string &gjfFileName);
    //! Create fake k.p atoms.
    void genKpAtoms(double Lx, double Ly, double ax, double ay, 
                    const vector<Atom>& periodicTable);

    // access functions
    Atom        AtomAt(uint i) const;       // Get one atom at i
    string      Symbol(int i) const { return convertIndexToSym(mAtomId(i)); } ;
    double      X(uint i) const { return mXyz(i, spacevec::X); };
    double      Y(uint i) const { return mXyz(i, spacevec::Y); };
    double      Z(uint i) const { return mXyz(i, spacevec::Z); };
    vec         X() const { return mXyz.col(spacevec::X); };
    vec         Y() const { return mXyz.col(spacevec::Y); };
    vec         Z() const { return mXyz.col(spacevec::Z); };
    int         NumOfAtoms() const { return mNumAtoms; };
    int         NumOfOrbitals() const { return mNumOrbitals; };
    int         NumOfElectrons() const { return mNumElectrons; };
    const lvec& LatticeVector() const { return mLatticeVector; };
    void        LatticeVector(const lvec& a) { this->mLatticeVector = a; };
    void        PeriodicTable(const vector<Atom> &periodicTable);
    double      xmin() {return min(mXyz.col(spacevec::X)); };
    double      xmax() {return max(mXyz.col(spacevec::X)); };
    double      ymin() {return min(mXyz.col(spacevec::Y)); };
    double      ymax() {return max(mXyz.col(spacevec::Y)); };
    double      zmin() {return min(mXyz.col(spacevec::Z)); };
    double      zmax() {return max(mXyz.col(spacevec::Z)); };


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

}

#endif	/* ATOMS_H */

