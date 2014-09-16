/* 
 * File:   AtomicStruct.h
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
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

#include "Lattice.h"

#include "maths/svec.h"
#include "maths/arma.hpp"
#include "grid/grid.hpp"
#include "string/stringutils.h"
#include "utils/Printable.hpp"
#include "utils/serialize.hpp"
#include "utils/vout.h"

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace qmicad{
namespace atoms{

using utils::Printable;
using utils::trim;
using utils::meshgrid; 
using utils::linspace;
using namespace utils::stds;
using namespace maths::armadillo;
using namespace maths::spvec;

using boost::shared_ptr;
using boost::make_shared;  
    
        
/**
 * The Atom structure describes properties of an atom as seen in 
 * the periodic table.
 */
struct Atom {
    uint    ia;   //!< Atomic number
    string  sym;  //!< Atomic symbol
    uint    ne;   //!< Number of electrons
    uint    no;   //!< Number of orbitals

    //! Default constructor.
    Atom() {};
    //! Detailed constructor taking all the properties.
    Atom(uint ia, string sym, uint ne, uint no){
        this->ia = ia;
        this->sym = sym;
        this->ne = ne;
        this->no = no;
    };
    
    //! Copy constructor.
    Atom(const Atom& orig): ia(orig.ia), sym(orig.sym), ne(orig.ne), no(orig.no){
    };
        
    //!< Compares everything.
    friend bool operator == (const Atom &lhs, const Atom &rhs){
        return (lhs.ia == rhs.ia && lhs.sym == rhs.sym &&
                lhs.ne == rhs.ne && lhs.no == rhs.no);
    }
    
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
 * Our periodic table.
 */
struct PeriodicTable{
    typedef map<uint, Atom>::const_iterator cpiter;
    typedef map<uint, Atom>::const_iterator piter;
    
    map<uint, Atom> elements;
    
    PeriodicTable(){
        init();
    };
    
    void init(){
        // Fill up the periodic table
        add(0,  "D",  1, 1); // For discretized Hamiltonian
        add(1,  "H",  1, 1);
        add(5,  "B",  1, 1);
        add(6,  "C",  1, 1);
        add(7,  "N",  1, 1);
        add(14, "Si", 1, 1);
        add(16, "S",  1, 1);
        add(32, "Ge", 1, 1);
        add(42, "Mo", 1, 1);
    }
    
    void add(uint ia, string sym, uint ne, uint no){
        elements[ia] = Atom(ia, sym, ne, no);
    }
    
    void add(const Atom &a){
        elements[a.ia] = a;
    }
    
    //!< Find atomic symbol.
    int find(const string& sym) const{
        for (cpiter it = elements.begin(); it != elements.end(); ++it){
            if (it->second.sym == sym){
                return it->second.ia;
            }
        }
        return -1; // not found
    }

    //!< Find atomic number
    int find(uint ia) const {
        cpiter it = elements.find(ia);
        if (it == elements.end()){ 
            return -1; // not found
        }else{
            return it->second.ia;
        }       
    }
    
    //!< Find atom.
    int find(const Atom& a) const{
        return find(a.ia);
    }

    int size() const{
        return elements.size();
    }
    
    Atom operator[](uint ia) const{
        return elements.find(ia)->second;
    }
    
    Atom operator[](const Atom & a) const{
        return elements.find(a.ia)->second;
    }

    Atom at(uint ia){
        return this->operator [](ia);
    }
    
    void update(const PeriodicTable pt){
        cpiter it;
        for(it = pt.elements.begin(); it != pt.elements.end(); ++it){
            // Copy and overwrite everything of new pt to this periodic table. 
            add(it->second);
        }
    }
    
    //! For MPI send/receive.
    friend class boost::serialization::access;

    //! For MPI send/receive.
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        ar & elements;
    }

};

typedef PeriodicTable ptable;

/** 
 * All atoms in the structure.
 */
class AtomicStruct: public Printable {
//typedefs
public:
    typedef shared_ptr<AtomicStruct> ptr;
// Fields
protected:
    ptable mpt;         //!< Our periodic table.
    int mNa;            //!< No of atoms.
    int mNo;            //!< No of orbitals
    int mNe;            //!< No of electrons

    // atoms and their coordinates
    icol mia;           //!< Atomic numbers for all atoms in the collection
    mat  mXyz;          //!< Atomic coordinates
    lvec mlv;           //!< Lattice vector

// Methods    
public:
    //!< Default constructor.
    AtomicStruct();
    //!< Copy constructor.
    AtomicStruct(const AtomicStruct& orig);
    //!< constructs from a GaussView GJF file.
    AtomicStruct(const string& gjfFileName ); 
    //!< constructs a periodic table.
    AtomicStruct(const ptable &periodicTable);    
    //!< constructs from a GaussView GJF file and a periodic table.
    AtomicStruct(const string& gjfFileName, const ptable &periodicTable);
    //!< Constructs from a atomic coordinates and lattice vector.
    AtomicStruct(const icol& atomId, const mat& coordinate, const lvec& lv);
    //!< Constructs from a atomic coordinates and lattice vector.
    AtomicStruct(const icol& atomId, const mat& coordinate, const lvec& lv,
    const ptable& periodicTable);
    //!< Destructor.
    virtual ~AtomicStruct(){};
    friend void swap(AtomicStruct& first, AtomicStruct& second);
    // operators
    AtomicStruct  operator()(maths::armadillo::span s) const;// Get a sub cell 
    AtomicStruct  operator()(const ucol& index) const;       // Get a sub cell
    AtomicStruct  operator()(uint i) const;                  // Get one atom cell
    
    AtomicStruct& operator= (AtomicStruct rhs);              // assignment
    AtomicStruct& operator+= (const AtomicStruct& atoms);    // concatenation
    AtomicStruct& operator-= (const lcoord& latticeCoord);   // coordinate shifting
    AtomicStruct& operator-= (const svec& positionVect);     // coordinate shifting
    AtomicStruct& operator+= (const lcoord& latticeCoord);   // coordinate shifting
    AtomicStruct& operator+= (const svec& positionVect);     // coordinate shifting

    // non-mamber operators
    friend AtomicStruct operator- (AtomicStruct atm, const lcoord& latticeCoord);     // atm2 = atm1 - lc;
    friend AtomicStruct operator- (AtomicStruct atm, const svec& positionVect);       // atm2 = atm1 - r;
    friend AtomicStruct operator+ (AtomicStruct atm, const lcoord& latticeCoord);     // atm2 = atm1 + lc;
    friend AtomicStruct operator+ (AtomicStruct atm, const svec& positionVect);       // atm2 = atm1 + r;
    friend AtomicStruct operator+ (AtomicStruct atmi, const AtomicStruct& atmj);             // concatanation

    // utilities
    //!< Import atoms from Gaussview file.
    void importGjf(const string &gjfFileName);
    //!< Export atoms to Gaussview file.
    void exportGjf(const string &gjfFileName);
    //!< Generate atoms in a rectangular lattice.
    void genRectLattAtoms(uint nl, uint nw, double ax, double ay, 
                    const ptable& periodicTable);
    void genSimpleCubicStruct(const Atom &atom, double a, uint nl, uint nw = 1, uint nh = 1);
    //!< Generates simple cubic structures.
    void genSimpleCubicStruct(const Atom &atom, double a, double l, double w = 0, double h = 0);
    //!< String representation of Atomic structure.
    string toString() const;

    // access functions
    Atom        AtomAt(uint i) const;       // Get one atom at i
    string      Symbol(uint i) const { return mpt[mia(i)].sym; } ;
    double      X(uint i) const { return mXyz(i, coord::X); };
    double      Y(uint i) const { return mXyz(i, coord::Y); };
    double      Z(uint i) const { return mXyz(i, coord::Z); };
    vec         X() const { return mXyz.col(coord::X); };
    vec         Y() const { return mXyz.col(coord::Y); };
    vec         Z() const { return mXyz.col(coord::Z); };
    mat         XYZ() const { return mXyz; };
    int         NumOfAtoms() const { return mNa; };
    int         NumOfOrbitals() const { return mNo; };
    int         NumOfElectrons() const { return mNe; };
    lvec        LatticeVector() const { return mlv; };
    void        LatticeVector(const lvec& a) { this->mlv = a; };
    void        PeriodicTable(const ptable &periodicTable);
    double      xmin() {return min(mXyz.col(coord::X)); };
    double      xmax() {return max(mXyz.col(coord::X)); };
    double      xl(){ return abs(xmax() - xmin()); };
    double      ymin() {return min(mXyz.col(coord::Y)); };
    double      ymax() {return max(mXyz.col(coord::Y)); };
    double      yl(){ return abs(ymax() - ymin()); };
    double      zmin() {return min(mXyz.col(coord::Z)); };
    double      zmax() {return max(mXyz.col(coord::Z)); };
    double      zl(){ return abs(zmax() - zmin()); };
    
    AtomicStruct span(uint start, uint end) const;


protected:
    void        init();
    int         computeNumOfOrbitals();
    int         computeNumOfElectrons();

    //!< Serialization.
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version){
        Printable::serialize(ar, version);
        ar & mpt;
        ar & mNa;
        ar & mNo;
        ar & mNe;
        ar & mia;
        ar & mXyz;
        ar & mlv;
    }
    
};


}
}

#endif	/* ATOMS_H */

