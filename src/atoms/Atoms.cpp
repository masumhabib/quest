/* 
 * File:   Atoms.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 * 
 * Created on April 5, 2013, 1:27 PM
 * 
 * Description: The class Atoms encapsulates everything about the atoms present in
 * the device.
 * 
 */

#include "Atoms.h"

namespace qmicad{
/* Default constructor */
Atoms::Atoms() {
    initPeriodicTable();
    init();
}

// constructs from a GaussView gjf file and a periodic table
Atoms::Atoms(const string& gjfFileName, const vector<Atom>& PeriodicTable):
mPeriodicTable(PeriodicTable){    
    init();
    importGjf(gjfFileName);
}

/* Construct from a GaussView gjf file */
Atoms::Atoms(const string& gjfFileName) {
    initPeriodicTable();
    init();
    
    importGjf(gjfFileName);
}

/* Copy constructor */
Atoms::Atoms(const Atoms& orig):
mLatticeVector(orig.mLatticeVector),
mAtomId(orig.mAtomId),
mPeriodicTable(orig.mPeriodicTable),
mXyz(orig.mXyz)
{
    mNumAtoms = orig.mNumAtoms;
    mNumOrbitals = orig.mNumOrbitals;
    mNumElectrons = orig.mNumElectrons;
}

/* Construct from coordinates */
Atoms::Atoms(const icol& atomId, const mat& coordinate, const lvec& lv)
{
    if (atomId.n_rows != coordinate.n_rows){
        throw invalid_argument("In Atoms::Atoms(const icol& atomId, const mat& "
                "coordinate, const lvec& lv) number of rows of atomId and "
                "coordinate does not match");
    }

    initPeriodicTable();
    
    mLatticeVector = lv;
    mAtomId = atomId;
    mXyz = coordinate;

    mNumAtoms = mAtomId.n_rows;
    mNumOrbitals = computeNumOfOrbitals();
    mNumElectrons = computeNumOfElectrons();    
}

/* Construct from coordinates and a periodic table*/
Atoms::Atoms(const icol& atomId, const mat& coordinate, const lvec& lv,
const vector<Atom>& PeriodicTable):mPeriodicTable(PeriodicTable){
    if (atomId.n_rows != coordinate.n_rows){
        throw invalid_argument("In Atoms::Atoms(const icol& atomId, const mat& "
                "coordinate, const lvec& lv) number of rows of atomId and "
                "coordinate does not match");
    }
    
    mLatticeVector = lv;
    mAtomId = atomId;
    mXyz = coordinate;

    mNumAtoms = mAtomId.n_rows;
    mNumOrbitals = computeNumOfOrbitals();
    mNumElectrons = computeNumOfElectrons();    
}

Atoms::Atoms(double Lx, double Ly, double ax, double ay, const vector<Atom>& 
    periodicTable)
{
    init();
    genKpAtoms(Lx, Ly, ax, ay, periodicTable);
}

/* our periodic table */
void Atoms::initPeriodicTable(){
    // Fill up the periodic table
    mPeriodicTable.push_back(Atom(0,  "D",  1, 1)); // For discretized Hamiltonian
    mPeriodicTable.push_back(Atom(1,  "H",  1, 1));
    mPeriodicTable.push_back(Atom(5,  "B",  1, 1));
    mPeriodicTable.push_back(Atom(6,  "C",  1, 1));
    mPeriodicTable.push_back(Atom(7,  "N",  1, 1));
    mPeriodicTable.push_back(Atom(14, "Si", 1, 1));
    mPeriodicTable.push_back(Atom(16, "S",  1, 1));
    mPeriodicTable.push_back(Atom(32, "Ge", 1, 1));
    mPeriodicTable.push_back(Atom(42, "Mo", 1, 1));

}

/* initializer */
void Atoms::init(){
    mNumAtoms = 0;
    mNumOrbitals = 0;
    mNumElectrons = 0;
    mAtomId.clear();
    mXyz.clear();
    mLatticeVector.zeros();
}

/* Assignment operator: 
 * c = b; creates a copy of b and assigns it to c*/
Atoms& Atoms::operator= (Atoms rhs){
    swap(*this, rhs);
    
    return *this;
}

/* Swaps two Atoms objects */
void swap(Atoms& first, Atoms& second){
    using std::swap;

    swap(first.mLatticeVector, second.mLatticeVector);
    swap(first.mAtomId, second.mAtomId);
    swap(first.mNumAtoms, second.mNumAtoms);
    swap(first.mNumElectrons, second.mNumElectrons);
    swap(first.mNumOrbitals, second.mNumOrbitals);
    swap(first.mPeriodicTable, second.mPeriodicTable);
    
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
Atoms operator+ (Atoms atmi, const Atoms& atmj){
    atmi += atmj;
    return atmi;
}

/* Concatenation: atmi += atmj */
Atoms& Atoms::operator+= (const Atoms& atj){

    // update the periodic table
    for(int it = 0; it != atj.mPeriodicTable.size(); ++it){
        // Update our table if we do not already have  
        // the atoms in atj;
        if (convertSymToIndex(atj.mPeriodicTable[it].sym) == -1){ // not found in our database
            mPeriodicTable.push_back(atj.mPeriodicTable[it]); // insert it to our database
        }
    }
    
    // concatenate the atom id's and coordinates
    mAtomId.insert_rows(mNumAtoms,atj.mAtomId);
    mXyz.insert_rows(mNumAtoms,atj.mXyz);
    
    // add the lattice vectors
    mLatticeVector += atj.mLatticeVector;
    
    // add the number of atoms, orbitals and electrons
    mNumAtoms += atj.mNumAtoms;
    mNumOrbitals += atj.mNumOrbitals;
    mNumElectrons += atj.mNumElectrons;
    
    return *this;
}

/* Atoms - Lattice coordinate*/
Atoms operator- (Atoms atm, const lcoord& lc){
    atm -= lc;
    return atm;
};

/* Atoms -= Lattice coordinate*/
Atoms& Atoms::operator-= (const lcoord& lc){
    
    *this -= mLatticeVector.a1*lc.n1;
    *this -= mLatticeVector.a2*lc.n2;
    *this -= mLatticeVector.a3*lc.n3;
    
    return *this;
}

/* Atoms = Atoms - position vector */
Atoms operator- (Atoms atm, const svec& r){
    atm -= r;
    return atm;
};

/* Atoms -= position vector */
Atoms& Atoms::operator-= (const svec& rvec){
    
    mXyz.col(spacevec::X) -= rvec(spacevec::X);
    mXyz.col(spacevec::Y) -= rvec(spacevec::Y);
    mXyz.col(spacevec::Z) -= rvec(spacevec::Z);
    
    return *this;
}

/* Atoms + Lattice coordinate*/
Atoms operator+ (Atoms atm, const lcoord& lc){
    atm += lc;
    return atm;
};

/* Atoms += Lattice coordinate*/
Atoms& Atoms::operator+= (const lcoord& lc){
    
    *this += mLatticeVector.a1*lc.n1;
    *this += mLatticeVector.a2*lc.n2;
    *this += mLatticeVector.a3*lc.n3;
    
    return *this;
}

/* Atoms = Atoms + position vector */
Atoms operator+ (Atoms atm, const svec& r){
    atm += r;
    return atm;
};

/* Atoms += position vector */
Atoms& Atoms::operator+= (const svec& rvec){
    
    mXyz.col(spacevec::X) += rvec(spacevec::X);
    mXyz.col(spacevec::Y) += rvec(spacevec::Y);
    mXyz.col(spacevec::Z) += rvec(spacevec::Z);
    
    return *this;
}

/*
 * Callable operators to extract a subset of atoms whose atom indices are stored
 * in column vector.
 * @FIXME: these operators does not guarantee correct lattice vector.
 */
Atoms Atoms::operator ()(const ucol& index) const{
    
    //get only the atoms we are interested in
    icol atomId = mAtomId.elem(index);   
    ucol cols;
    cols << spacevec::X << spacevec::Y << spacevec::Z;
    mat coordinate = mXyz(index,cols);        

    return Atoms(atomId, coordinate, mLatticeVector, mPeriodicTable);
}

/*
 * Callable operators to extract a subset of atoms defined by a span.
 * @FIXME: these operators does not guarantee correct lattice vector.
 */
Atoms Atoms::operator ()(span s) const{
    
    //get only the atoms we are interested in
    icol atomId = mAtomId(s);   
    mat coordinate = mXyz(s,span::all); 

    return Atoms(atomId, coordinate, mLatticeVector, mPeriodicTable);
}

/*
 * Callable operators to extract a one atom with index i.
 * @FIXME: these operators does not guarantee correct lattice vector.
 */

Atoms Atoms::operator ()(uint i) const{
    
    ucol index(1);
    index(0) = i;
    return (*this)(index);
}

Atom Atoms::AtomAt(uint i) const{
    return mPeriodicTable[mAtomId[i]];
}

/* Dump the data to the stream */
ostream& operator << (ostream & out, const Atoms &b)
{
    
    // write the file header
    out << "# hf/3-21g" << endl; // header command
    out << endl;
    out << "GRAPHENE" << endl;   // tittle
    out << endl;
    out << "0 1" << endl;        // charge spin
    
    // write atomic coordinates
    out.precision(4);
    out << std::fixed;
    for(int i=0; i < b.NumOfAtoms(); i++){
        out.width(3);
        out << b.Symbol(i);
        out.width(10); 
        out << std::fixed << b.X(i);
        out.width(10);
        out << b.Y(i);
        out.width(10);
        out << b.Z(i) << endl;
    }
    
    out << endl;
    
	return out;
}

void Atoms::importGjf(const string& gjfFileName){
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

            ind = convertSymToIndex(sym);
            if (ind > -1){
                
                mAtomId.resize(mNumAtoms+1);
                mAtomId(mNumAtoms) =  ind;
                
                row r(3);
                r << x << y << z;
                mXyz.insert_rows(mXyz.n_rows,r);
                
                mNumOrbitals += mPeriodicTable[ind].no;
                mNumElectrons += mPeriodicTable[ind].ne;
                mNumAtoms++;
            }

        // If we find the line that contains "0 1" which is 
        // the spin and charge specification, then we are in the 
        // coordinate section of the gjf file.
        }else if(inHeader && isdigit(line[0])){ 
            inHeader = false;
        }        
    }
    
    if (mNumAtoms == 0){
        mAtomId.reset();
        mXyz.reset();
        
        throw runtime_error(" No atoms found in " + gjfFileName + ".");;
    }
}

void Atoms::exportGjf(const string& gjfFileName){
    
    ofstream gjf(gjfFileName.c_str());
    
    if (!gjf.is_open()){
        throw runtime_error(" Failed to open file " + gjfFileName + ".");
    }
        
    gjf << *this;
    gjf.close();
}

void Atoms::genKpAtoms(double Lx, double Ly, double ax, double ay, 
        const vector<Atom>& periodicTable){
    MatGrid xy(-Lx/2, Lx/2, ax, -Ly/2, Ly/2, ay);
    
    // calculate x, y and z coordinates of the atoms
    mNumAtoms = xy.Nx()*xy.Ny();            // total number of atoms
    mXyz.set_size(mNumAtoms, 3);            // xyz coordinate of atoms
    mXyz.col(spacevec::X) = xy.X();
    mXyz.col(spacevec::Y) = xy.Y();
    mXyz.col(spacevec::Z).zeros();
    
    // prepare atomId list containing atomic number of a fake atom 'D'.
    mAtomId.set_size(mNumAtoms);
    PeriodicTable(periodicTable);           // Copy the periodic table.
    mAtomId.fill(mPeriodicTable[0].ia);
    
    // lattice vector
    mLatticeVector.a1(spacevec::X) = xy.maxx()-xy.minx() + ax;
    mLatticeVector.a2(spacevec::Y) = xy.maxy()-xy.miny() + ay;  
    
    // calculate number of orbitals and electrons.
    mNumOrbitals = computeNumOfOrbitals();
    mNumElectrons = computeNumOfElectrons();    
 
}

int Atoms::convertSymToIndex(const string& sym) const{
    for (int it = 0; it < mPeriodicTable.size(); ++it){
        if (mPeriodicTable[it].sym == sym){
            return it;
        }
    }
    
    return -1; // not found
}

string Atoms::convertIndexToSym(int ind) const {
    return mPeriodicTable[ind].sym;
}

int Atoms::computeNumOfOrbitals(){
    int numOrbitals = 0;
    for(int ia = 0; ia < mNumAtoms; ++ia){
        numOrbitals += mPeriodicTable[mAtomId[ia]].no; 
    }
    
    return numOrbitals;
}

int Atoms::computeNumOfElectrons(){
    int numElectrons = 0;
    for(int ia = 0; ia < mNumAtoms; ++ia){
        numElectrons += mPeriodicTable[mAtomId[ia]].ne;       
    }
    
    return numElectrons;
}

void Atoms::PeriodicTable(const vector<Atom>& periodicTable){
    mPeriodicTable = periodicTable; 
    mNumOrbitals = computeNumOfOrbitals();
    mNumElectrons = computeNumOfElectrons();
}

}


