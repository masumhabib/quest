/* 
 * File:   AtomicStruct.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 5, 2013, 1:27 PM
 * 
 * Description: The class Atoms encapsulates everything about the atoms present in
 * the device structure.
 * 
 */

#include "atoms/AtomicStruct.h"

namespace qmicad{
namespace atoms{

/* Default constructor */
AtomicStruct::AtomicStruct():mpt() {
    init();
}

/* Copy constructor */
AtomicStruct::AtomicStruct(const AtomicStruct& orig):
mpt(orig.mpt),
mia(orig.mia),
mXyz(orig.mXyz),
mlv(orig.mlv)
{
    mNa = orig.mNa;
    mNo = orig.mNo;
    mNe = orig.mNe;
}

/* Construct from a GaussView gjf file */
AtomicStruct::AtomicStruct(const string& gjfFileName):mpt() {
    init();
    
    importGjf(gjfFileName);
}

// constructs from a GaussView gjf file and a periodic table
AtomicStruct::AtomicStruct(const string& gjfFileName, const ptable &periodicTable):
mpt(periodicTable){    
    init();
    importGjf(gjfFileName);
}

// constructs from a periodic table
AtomicStruct::AtomicStruct(const ptable &periodicTable):
mpt(periodicTable){    
    init();
}

/* Construct from coordinates */
AtomicStruct::AtomicStruct(const icol& atomId, const mat& coordinate, const lvec& lv):
mpt()
{
    if (atomId.n_rows != coordinate.n_rows){
        throw invalid_argument("In Atoms::Atoms(const icol& atomId, const mat& "
                "coordinate, const lvec& lv) number of rows of atomId and "
                "coordinate does not match");
    }
    
    mlv = lv;
    mia = atomId;
    mXyz = coordinate;

    mNa = mia.n_rows;
    mNo = computeNumOfOrbitals();
    mNe = computeNumOfElectrons();    
}

/* Construct from coordinates and a periodic table*/
AtomicStruct::AtomicStruct(const icol& atomId, const mat& coordinate, const lvec& lv,
const ptable& periodicTable):mpt(periodicTable){
    if (atomId.n_rows != coordinate.n_rows){
        throw invalid_argument("In Atoms::Atoms(const icol& atomId, const mat& "
                "coordinate, const lvec& lv) number of rows of atomId and "
                "coordinate does not match");
    }
    
    mlv = lv;
    mia = atomId;
    mXyz = coordinate;

    mNa = mia.n_rows;
    mNo = computeNumOfOrbitals();
    mNe = computeNumOfElectrons();    
}


/* initializer */
void AtomicStruct::init(){
    mNa = 0;
    mNo = 0;
    mNe = 0;
    mia.clear();
    mXyz.clear();
    mlv.zeros();
    mlv.Prefix(" ");
}

/* Assignment operator: 
 * c = b; creates a copy of b and assigns it to c*/
AtomicStruct& AtomicStruct::operator= (AtomicStruct rhs){
    swap(*this, rhs);
    
    return *this;
}

/* Swaps two Atoms objects */
void swap(AtomicStruct& first, AtomicStruct& second){
    using std::swap;

    swap(first.mlv, second.mlv);
    swap(first.mia, second.mia);
    swap(first.mNa, second.mNa);
    swap(first.mNe, second.mNe);
    swap(first.mNo, second.mNo);
    swap(first.mpt, second.mpt);
    
    /* The swap function in armadillo probably has a bug
     * that prevents swapping a zero sized matrix
     */
    if (first.mXyz.empty() || second.mXyz.empty()){
        std::swap(first.mXyz, second.mXyz);
    }else{
        swap(first.mXyz, second.mXyz);
    }
    
}

/* Concatenation: atmi + atmj */
AtomicStruct operator+ (AtomicStruct atmi, const AtomicStruct& atmj){
    atmi += atmj;
    return atmi;
}

/* Concatenation: atmi += atmj */
AtomicStruct& AtomicStruct::operator+= (const AtomicStruct& atj){
    // update our periodic table
    mpt.update(atj.mpt);
    
    // concatenate the atom id's and coordinates
    mia.insert_rows(mNa,atj.mia);
    mXyz.insert_rows(mNa,atj.mXyz);
    
    // add the lattice vectors
    mlv += atj.mlv;
    
    // add the number of atoms, orbitals and electrons
    mNa += atj.mNa;
    mNo += atj.mNo;
    mNe += atj.mNe;
    
    return *this;
}

/* Atoms - Lattice coordinate*/
AtomicStruct operator- (AtomicStruct atm, const lcoord& lc){
    atm -= lc;
    return atm;
};

/* Atoms -= Lattice coordinate*/
AtomicStruct& AtomicStruct::operator-= (const lcoord& lc){
    
    *this -= mlv.a1*lc.n1;
    *this -= mlv.a2*lc.n2;
    *this -= mlv.a3*lc.n3;
    
    return *this;
}

/* Atoms = Atoms - position vector */
AtomicStruct operator- (AtomicStruct atm, const svec& r){
    atm -= r;
    return atm;
};

/* Atoms -= position vector */
AtomicStruct& AtomicStruct::operator-= (const svec& rvec){
    
    mXyz.col(coord::X) -= rvec(coord::X);
    mXyz.col(coord::Y) -= rvec(coord::Y);
    mXyz.col(coord::Z) -= rvec(coord::Z);
    
    return *this;
}

/* Atoms + Lattice coordinate*/
AtomicStruct operator+ (AtomicStruct atm, const lcoord& lc){
    atm += lc;
    return atm;
};

/* Atoms += Lattice coordinate*/
AtomicStruct& AtomicStruct::operator+= (const lcoord& lc){
    
    *this += mlv.a1*lc.n1;
    *this += mlv.a2*lc.n2;
    *this += mlv.a3*lc.n3;
    
    return *this;
}

/* Atoms = Atoms + position vector */
AtomicStruct operator+ (AtomicStruct atm, const svec& r){
    atm += r;
    return atm;
};

/* Atoms += position vector */
AtomicStruct& AtomicStruct::operator+= (const svec& rvec){
    
    mXyz.col(coord::X) += rvec(coord::X);
    mXyz.col(coord::Y) += rvec(coord::Y);
    mXyz.col(coord::Z) += rvec(coord::Z);
    
    return *this;
}

/*
 * Callable operators to extract a subset of atoms whose atom indices are stored
 * in column vector.
 * @FIXME: these operators does not guarantee correct lattice vector.
 */
AtomicStruct AtomicStruct::operator ()(const ucol& index) const{
    
    //get only the atoms we are interested in
    icol atomId = mia.elem(index);   
    ucol cols;
    cols << coord::X << coord::Y << coord::Z;
    mat coordinate = mXyz(index,cols);        

    return AtomicStruct(atomId, coordinate, mlv, mpt);
}

/*
 * Callable operators to extract a subset of atoms defined by a span.
 * @FIXME: these operators does not guarantee correct lattice vector.
 */
AtomicStruct AtomicStruct::operator ()(maths::armadillo::span s) const{
    
    //get only the atoms we are interested in
    icol atomId = mia(s);   
    mat coordinate = mXyz(s,span::all); 

    return AtomicStruct(atomId, coordinate, mlv, mpt);
}

/*
 * Callable operators to extract a one atom with index i.
 * @FIXME: these operators does not guarantee correct lattice vector.
 */
AtomicStruct AtomicStruct::operator ()(uint i) const{
    
    ucol index(1);
    index(0) = i;
    return (*this)(index);
}

Atom AtomicStruct::AtomAt(uint i) const{
    return mpt[mia[i]];
}

string AtomicStruct::toString() const{
    
    stringstream out;
    // write atomic coordinates
    out << " Atomic Structure: " << endl;
    out.precision(4);
    out << std::fixed;
    for(int i=0; i < NumOfAtoms(); i++){
        out.width(3);
        out << Symbol(i);
        out.width(15); 
        out << std::fixed << X(i);
        out.width(15);
        out << Y(i);
        out.width(15);
        out << Z(i) << endl;
    }
    out << endl;
    
    // Lattice vectors.
    out << " Lattice Vector: " << endl;
    out << mlv << endl;
    
    // other information
    out << " Total atoms: " << mNa << endl;
    out << " Total orbitals: " << mNo << endl;
    out << " Total electrons: " << mNe << endl;
    
	return out.str();
}


void AtomicStruct::importGjf(const string& gjfFileName){
    string line;
    ifstream gjf (gjfFileName.c_str());        
    
    if(!gjf.is_open()){
        throw runtime_error("Failed to open file " + gjfFileName + ".");
    }
    
    // reset values
    init();
    
    bool inHeader = true;
    while(gjf.good()){
        getline(gjf, line);

        line = trim(line); 
        // Skip the comment, the command and the blank lines
        // and get the coordinates
        if(!inHeader && line.length() > 0 && line[0] != '#'
           && line[0] != '%' && line[0] != '\n'){
            
            stringstream ssline (line);
            string sym;
            double x, y, z;
            int ind;
            
            ssline >> sym >> x >> y >> z;

            ind = mpt.find(sym);
            if (ind > -1){
                
                mia.resize(mNa+1);
                mia(mNa) =  ind;
                
                row r(3);
                r << x << y << z;
                mXyz.insert_rows(mXyz.n_rows,r);
                
                mNo += mpt[ind].no;
                mNe += mpt[ind].ne;
                mNa++;
            }

        // If we find the line that contains "0 1" which is 
        // the spin and charge specification, then we are in the 
        // coordinate section of the gjf file.
        }else if(inHeader && isdigit(line[0])){ 
            inHeader = false;
        }        
    }
    
    if (mNa == 0){
        mia.reset();
        mXyz.reset();
        
        throw runtime_error(" No atoms found in " + gjfFileName + ".");;
    }
}

void AtomicStruct::exportGjf(const string& gjfFileName){
    
    ofstream gjf(gjfFileName.c_str());
    
    if (!gjf.is_open()){
        throw runtime_error(" Failed to open file " + gjfFileName + ".");
    }
    
    // write the file header
    gjf << "# hf/3-21g" << endl; // header command
    gjf << endl;
    gjf << "Generated by QMICAD" << endl;   // tittle
    gjf << endl;
    gjf << "0 1" << endl;        // charge spin
    
    // write atomic coordinates
    gjf.precision(4);
    gjf << std::fixed;
    for(int i=0; i < NumOfAtoms(); i++){
        gjf.width(3);
        gjf << Symbol(i);
        gjf.width(15); 
        gjf << std::fixed << X(i);
        gjf.width(15);
        gjf << Y(i);
        gjf.width(15);
        gjf << Z(i) << endl;
    }
    gjf << endl;
    
    gjf.close();
}

void AtomicStruct::genRectLattAtoms(uint nl, uint nw, double ax, double ay, 
        const ptable &periodicTable){
    
    vout << endl << " Deprecation Warning: AtomicStruct::genRectLattAtoms() will be removed from future release.";
    vout << " Use AtomicStruct::genSimpleCubicStruct() instead." << endl;
    
    double w = (nw-1)*ax;                  // width (in A)
    double l = (nl-1)*ay;                  // length  (in A)

    // create meshgrid of k.p atoms    
    col X, Y;
    meshgrid(X, Y, -l/2, l/2, ax, -w/2, w/2, ay);
    
    // calculate x, y and z coordinates of the atoms
    mNa = X.n_rows;            // total number of atoms
    mXyz.set_size(mNa, 3);            // xyz coordinate of atoms
    mXyz.col(coord::X) = X;
    mXyz.col(coord::Y) = Y;
    mXyz.col(coord::Z).zeros();
    
    // prepare atomId list containing atomic number of a fake atom 'D'.
    mia.set_size(mNa);
    mpt = periodicTable;           // Copy the periodic table.
    mia.fill(mpt[0].ia);
    
    // lattice vector
    mlv.a1(coord::X) = max(X)-min(X) + ax;
    mlv.a2(coord::Y) = max(Y)-min(Y) + ay;  
    
    // calculate number of orbitals and electrons.
    mNo = computeNumOfOrbitals();
    mNe = computeNumOfElectrons();    
 
}

void AtomicStruct::genSimpleCubicStruct(const Atom &atom, double a, double l, double w, double h){
    // create grid along x, y and z axes.
    col X, Y, Z;
    X = linspace<double>(-l/2, l/2, a);
    Y = linspace<double>(-w/2, w/2, a);
    Z = linspace<double>(-h/2, h/2, a);
    
    // total number of atoms
    long nx = X.n_rows, ny = Y.n_rows, nz = Z.n_rows;
    mNa = nx*ny*nz;
        
    // prepare atomId list containing atomic number of atom.
    mpt.add(atom);
    mia.set_size(mNa);
    mia.fill(atom.ia);

    // create simple cubic lattice.
    mXyz.set_size(mNa, 3);            // xyz coordinate of atoms
    long ia = 0;
    for (long ix = 0; ix < nx; ++ ix){
        for (long iy = 0; iy < ny; ++iy){
            for(long iz = 0; iz < nz; ++iz){
                mXyz(ia, coord::X) = X(ix);
                mXyz(ia, coord::Y) = Y(iy);
                mXyz(ia, coord::Z) = Z(iz);                
                ++ia;
            }
        }
    }
    
    // update lattice vector
    mlv.a1(coord::X) = max(X)-min(X) + a;
    mlv.a2(coord::Y) = max(Y)-min(Y) + a;  
    mlv.a3(coord::Z) = max(Z)-min(Z) + a;  
    
    // calculate number of orbitals and electrons.
    mNo = computeNumOfOrbitals();
    mNe = computeNumOfElectrons();
}

void AtomicStruct::genGNR(const Atom &atom, double acc, double l, double w, double h)
{
    // No of primitive cell needed in x direction
    double nbx = l/(3*acc);
    
    if( ceil(nbx)*3*acc - acc   <=   l ){
        nbx = ceil(nbx);
    }else{
        nbx = floor(nbx);
    }
    // No of primitive cell needed in y direction
    double nby = w / ( sqrt(3)*acc );
    
    if(  ceil(nby) * sqrt(3) * acc  -  sqrt(3) * acc / 2     <=    w  ){
        nby = ceil(nby);
    }else{
        nby = floor(nby);
    }
    // generating primitive cell
    AtomicStruct basisStructForGNR = this->genGNRPrimitiveCell( atom, acc );
    
    AtomicStruct wholeGNR;
    wholeGNR.init();
    
    // Creating the GNR structure using the primitive cell and its lattice vector
    for (long ix = 0; ix < nbx; ++ix){
        for ( long iy=0; iy < nby; ++iy ){
            AtomicStruct tempBasisStruct = basisStructForGNR;
            
            tempBasisStruct += ix * basisStructForGNR.mlv.a1; // shifting basisStruct
            tempBasisStruct += iy * basisStructForGNR.mlv.a2; // shifting basisStruct

            wholeGNR += tempBasisStruct; //concatenation 
        }
    }
    
    // overriding the lattice vectors // TODO // AtomicStruct operator += overloading adds up the lattice vector 
    wholeGNR.mlv.a1 =  basisStructForGNR.mlv.a1;
    wholeGNR.mlv.a2 =  basisStructForGNR.mlv.a2;
    wholeGNR.mlv.a3 =  basisStructForGNR.mlv.a3;
    
    // Assigning wholeGNR to this. 
    *this = wholeGNR; 
}

AtomicStruct AtomicStruct::genGNRPrimitiveCell(const Atom &atom, double acc){
    //////  generating the Primitive Cell consisting 4 atom for GNR
    //////                O      O
    //////           O                O
    //////
    
    // Updating Periodic Table
    mpt.add(atom);
    
    double nXorigin = 0;
    double nYorigin = 0;
    
    mat mCoOrdinates = zeros<mat>(4,3);
    mCoOrdinates( 0, coord::X ) = nXorigin;
    mCoOrdinates( 0, coord::Y ) = nYorigin;
    mCoOrdinates( 1, coord::X ) = nXorigin    +   acc / 2;
    mCoOrdinates( 1, coord::Y ) = nYorigin    +   sqrt(3) * acc / 2;
    mCoOrdinates( 2, coord::X ) = nXorigin    +   3 * acc / 2;
    mCoOrdinates( 2, coord::Y ) = nYorigin    +   sqrt(3) * acc / 2;
    mCoOrdinates( 3, coord::X ) = nXorigin    +   2 * acc;
    mCoOrdinates( 3, coord::Y ) = nYorigin;
    
    lvec tempMlv;
    tempMlv.zeros();
    tempMlv.a1(coord::X) = 3*acc;
    tempMlv.a2(coord::Y) = sqrt(3) * acc;
    
    icol tempAtomId = zeros<icol>(4);
    tempAtomId.fill(atom.ia);
    
    return AtomicStruct( tempAtomId, mCoOrdinates, tempMlv );

}

void AtomicStruct::genSimpleCubicStruct(const Atom &atom, double a, uint nl, uint nw, uint nh){
    double w = (nw-1)*a;                  // width (in A)
    double l = (nl-1)*a;                  // length  (in A)
    double h = (nh-1)*a;                  // length  (in A)
    
    genSimpleCubicStruct(atom, a, l, w, h);
}

int AtomicStruct::computeNumOfOrbitals(){
    int numOrbitals = 0;
    for(int ia = 0; ia < mNa; ++ia){
        numOrbitals += mpt[mia[ia]].no; 
    }
    
    return numOrbitals;
}

int AtomicStruct::computeNumOfElectrons(){
    int numElectrons = 0;
    for(int ia = 0; ia < mNa; ++ia){
        numElectrons += mpt[mia[ia]].ne;       
    }
    
    return numElectrons;
}

void AtomicStruct::PeriodicTable(const ptable& periodicTable){
    mpt = periodicTable; 
    mNo = computeNumOfOrbitals();
    mNe = computeNumOfElectrons();
}

AtomicStruct AtomicStruct::span(uint start, uint end) const{
    return this->operator()(maths::armadillo::span(start, end));
}

}
}

