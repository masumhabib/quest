/* 
 * File:   genHam.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 6, 2013, 5:52 PM
 * 
 * Description: Tight binding calculation logic and data.
 * 
 */

#ifndef HAMILTONIAN_H
#define	HAMILTONIAN_H

#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <armadillo>

#include "../atoms/AtomicStruct.h"
#include "../utils/vout.h"

namespace qmicad{
namespace hamiltonian{

using boost::shared_ptr;
using boost::make_shared;
using namespace utils::stds;
using utils::Printable;
using namespace maths::armadillo;
using atoms::AtomicStruct;

struct HamParams: public Printable{
    // Parameters required for all Hamiltonin
    double dtol;         // distance tolerance    
    
    HamParams(const string &prefix = ""):Printable(" " + prefix){
        mTitle = "Hamiltonian paramters";
        // default parameters          
    }
    
    // Updates internal parameters. Call it after changing any of the 
    // public parameters.
    virtual void update(){}; 
        
protected:

};

template<class T>
class Hamiltonian{
public:
    Hamiltonian(): mOrth(true){
    };
    
    Hamiltonian(const HamParams& hp): mOrth(true){
        mhp = (new HamParams(hp));
    };

    //!< Sets size of the internal storage of H and S matrices for 
    //!< NEGF calculations.
    virtual void setSizeForNegf(uint nbi){
        setSize(nbi, 2);
    }

    //!< Sets size of the internal storage of H and S matrices for
    //!< band structure calculation.
    virtual void setSizeForBand(uint nbi){
        setSize(nbi, 1);
    }

   //!< Sets size of the internal storage of H and S matrices.
    virtual void setSize(uint nbi, uint nbj){
        setHsize(nbi, nbj);
        setSsize(nbi, nbj);
    }

    //!< Sets size of the internal storage of H matrix. This does not allocate 
    //!< memory. Memory is allocated when the matrices are calculated.
    virtual void setHsize(uint nbi, uint nbj){
        mH.set_size(nbi, nbj);
    }
    //!< Sets size of the internal storage of S matrix.
    virtual void setSsize(uint nbi, uint nbj){
        mS.set_size(nbi, nbj);
    }

    //!< Returns the diagonal block of H for NEGF.
    shared_ptr<T> H0(uint ib){
        return mH(ib, 0);
    } 

    //!< Returns the lower diagonal block of H for NEGF.
    shared_ptr<T> Hl(uint ib){
        return mH(ib, 1);
    } 
    
    //!< Returns block (ib,jb) of H.
    shared_ptr<T> H(uint ib, uint jb = 0){
        return mH(ib, jb);
    } 

    //!< Returns a diagonal block of overlap matrix.
    shared_ptr<T> S0(uint ib){
        return mS(ib, 0);
    } 

    //!< Returns a lower diagonal block of overlap matrix.
    shared_ptr<T> Sl(uint ib){
        return mS(ib, 1);
    } 
    
    //!< Returns block (ib,jb) of overlap matrix.
    shared_ptr<T> S(uint ib, uint jb = 0){
        return mS(ib, jb);
    } 
    

    //!< Generates the hamiltonaina and overlap matrices along the diagonal block
    //!< and store them in (ib, 0) location.
    virtual void genDiagBlock(const AtomicStruct &bi, const AtomicStruct &bj, 
        uint ib)
    {
        generate(bi, bj, ib, 0);
    }

    //!< Generates the hamiltonaina and overlap matrices along the lower diagonal block
    //!< and store them in (ib, 0) location.
    virtual void genLowDiagBlock(const AtomicStruct &bi, const AtomicStruct &bj, 
        uint ib)
    {
        generate(bi, bj, ib, 1);
    }
    
    //!< Generates the hamiltonaina and overlap matrices between neighbors
    //!< and store them in (ib, 0) location.
    virtual void genNearestNeigh(const AtomicStruct &bi, const AtomicStruct &bj, 
        uint ib)
    {
        generate(bi, bj, ib, 0);
    }

    
    //!< Generates the hamiltonaina and overlap matrices between atomc block
    //!< i and atomic block j and store them in (ib, jb) location.
    virtual void generate(const AtomicStruct &bi, const AtomicStruct &bj, 
                          uint ib, uint jb)
    {
        // In orthogonal basis, only S_i,i are non zero.
        if (mOrth && ib != jb){
                T smat;
                mH(ib, jb) = make_shared<T>();
                generate(*mH(ib, jb), smat, bi, bj);
        }else{
            mH(ib, jb) = make_shared<T>();
            mS(ib, jb) = make_shared<T>();
            generate(*mH(ib, jb), *mS(ib, jb), bi, bj);
        }
    }

protected:    
    //!< Generate Hamiltonian between two atoms.
    virtual T genTwoAtomHam(const AtomicStruct& atomi, 
                            const AtomicStruct& atomj){};
    //!< Generate Overlap matrix between two atoms.
    virtual T genTwoAtomOvl(const AtomicStruct& atomi, 
                            const AtomicStruct& atomj){};

    //!< Generates the hamiltonaina and overlap matrices between atomc block
    //!< i and atomic block j. 
    virtual void generate(T &hmat, T&smat, const AtomicStruct &bi, const AtomicStruct &bj){
        // Just for easy reference
        int nai = bi.NumOfAtoms();
        int naj = bj.NumOfAtoms();
        int noi = bi.NumOfOrbitals();
        int noj = bj.NumOfOrbitals();

        // Most of the matrix elements are zeros. So, we'll only change the 
        // non zero elements below.
        hmat = zeros<T>(noi, noj);
        smat = zeros<T>(noi, noj);
        
        // Lets find the neighbors. 
        int io = 0;
        for(int ia = 0; ia != nai; ++ia){
            AtomicStruct atomi = bi(ia);        // extract atom ia
            int ni = atomi.NumOfOrbitals();     // number of orbitals in atom ia

            int jo = 0;
            for(int ja = 0; ja != naj; ++ja){
                AtomicStruct atomj = bj(ja);        // extract atom ja    
                int nj = atomj.NumOfOrbitals();     // number of orbitals in atom ja
                // generate Hamiltonian matrix between orbitals of
                // atom i and atom j
                hmat(span(io,io+ni-1), span(jo,jo+nj-1)) = genTwoAtomHam(atomi, atomj);
                smat(span(io,io+ni-1), span(jo,jo+nj-1)) = genTwoAtomOvl(atomi, atomj);
                jo += nj;
            }

            io += ni;
        }    
    }   

protected:
    shared_ptr<HamParams>  mhp;       //!< Hamiltonian generation parameters.
    
    bool                   mOrth;     //!< Is this an orthogonal basis.
    field<shared_ptr<T> >  mH;        //!< Block matrices of H: H_i,j.
    field<shared_ptr<T> >  mS;        //!< Block matrices of S: S_i,j.
};

typedef Hamiltonian<cxmat>  cxham;
typedef Hamiltonian<mat>    ham;

}
}
#endif	/* HAMILTONIAN_H */

