/* 
 * File:   genHam.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 6, 2013, 5:52 PM
 * 
 * Description: Tight binding calculation logic and data.
 * 
 */

#ifndef GENHAM_H
#define	GENHAM_H

#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <armadillo>

#include "../atoms/Atoms.h"

template<class T>class HamGenerator;

template<class T>
T genHam(const Atoms &ati, const Atoms &atj, const HamGenerator<T>& hg){
    // Just for easy reference
    int nai = ati.NumOfAtoms();
    int naj = atj.NumOfAtoms();
    int noi = ati.NumOfOrbitals();
    int noj = atj.NumOfOrbitals();

    // Most of the matrix elements are zeros. So, we'll only change the 
    // non zero elements below.
    T hmat = zeros<T>(noi, noj);
    
    // Lets find the neighbors. 
    int io = 0;
    for(int ia = 0; ia != nai; ++ia){
        Atoms atomi = ati(ia);              // extract atom ia
        int ni = atomi.NumOfOrbitals();     // number of orbitals in atom ia
        
        int jo = 0;
        for(int ja = 0; ja != naj; ++ja){
            Atoms atomj = atj(ja);              // extract atom ja    
            int nj = atomj.NumOfOrbitals();     // number of orbitals in atom ja
            // generate Hamiltonian matrix between orbitals of
            // atom i and atom j
            hmat(span(io,io+ni-1), span(jo,jo+nj-1)) = hg(atomi, atomj);
            jo += nj;
        }
        
        io += ni;
    }    
    return hmat;
}

template<class T>
class HamGenerator{
public:
    HamGenerator(){};
protected:
    virtual T operator()(const Atoms& atomi, const Atoms& atomj)const{};

public:
friend T genHam<T>(const Atoms &ati, const Atoms &atj, const HamGenerator<T>& pg);
    
};


#endif	/* GENHAM_H */

