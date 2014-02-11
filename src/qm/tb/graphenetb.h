/* 
 * File:   graphenetb.h
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on January 25, 2014, 1:04 AM
 * 
 * Tight binding parameters for graphene.
 */

#ifndef GRAPHENETB_HPP
#define	GRAPHENETB_HPP

#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>
#include <armadillo>

#include "../../maths/arma.hpp"
#include "../../atoms/AtomicStruct.h"
#include "../hamiltonian.hpp"

namespace qmicad{
using boost::static_pointer_cast;
using namespace maths::armadillo;
    
struct GrapheneTbParams: public HamParams{
    // Tight binding parameters
    float dtol;          // distance tolerance

    float ec;            // onside energy
    float di0;           // in-plane C-C bond length
    float ti0;           // in-plane C-C hopping
    float do0;           // out-of-plane C-C distance
    float to0;           // out-of-plane C-C hopping
    float doX;           // out-of-plane C-C neighbor cut-off distance in
                         // units of 1 C-C bond length
    // Inter layer all nearest neighbors
    // See PRL 109, 236604 (2012)
    float lmdz;         // lambda_z
    float lmdxy;        // lambda_xy
    float alpha;        // alpha
    
    // Constructor
    GrapheneTbParams(string prefix = ""):HamParams(prefix){
        mTitle = "Graphene";
    }
        
    string toString() const { 
        stringstream ss;
        ss << mPrefix << " dtol = " << dtol << endl;
        ss << mPrefix << " ec = " << ec << endl;
        ss << mPrefix << " di0 = " << di0 << ", ti0 = " << ti0 << endl;
        ss << mPrefix << " do0 = " << do0 << ", to0 = " << to0 
                      << ", doX = " << doX << endl;
        ss << mPrefix << " lmdz = " << lmdz << ", lmdxy = " << lmdxy 
                      << ", alpha = " << alpha;

        return ss.str(); 
    };   

};
/** 
 * Tight binding Hamiltonian and m for graphene.
 */
class GrapheneTbHam: public Hamiltonian<mat>{
// Methods
public:
    GrapheneTbHam(const GrapheneTbParams& p){
        mhp = shared_ptr<HamParams> (new GrapheneTbParams(p));
    };
    virtual ~GrapheneTbHam(){};

    //!< Generate Hamiltonian between two atoms.
    virtual mat genTwoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj);    
    //!< Generate overlap matrix between two atoms.
    virtual mat genTwoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj);    
};

}
#endif	/* GRAPHENETB_HPP */

