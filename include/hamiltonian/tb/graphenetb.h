/* 
 * File:   graphenetb.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 25, 2014, 1:04 AM
 * 
 * Tight binding parameters for graphene.
 */

#ifndef GRAPHENETB_HPP
#define	GRAPHENETB_HPP

#include <boost/serialization/string.hpp>
#include <boost/serialization/access.hpp>

#include "maths/arma.hpp"
#include "atoms/AtomicStruct.h"
#include "hamiltonian/hamiltonian.hpp"

namespace qmicad{
namespace hamiltonian{

using boost::static_pointer_cast;
using namespace maths::armadillo;

class GrapheneTbHam;

struct GrapheneTbParams: public HamParams{
    friend class GrapheneTbHam;
    // Tight binding parameters
    double ec;            // onside energy
    double di0;           // in-plane C-C bond length
    double ti0;           // in-plane C-C hopping
    double do0;           // out-of-plane C-C distance
    double to0;           // out-of-plane C-C hopping
    double doX;           // out-of-plane C-C neighbor cut-off distance in
                          // units of 1 C-C bond length
    // Inter layer all nearest neighbors
    // See PRL 109, 236604 (2012)
    double lmdz;         // lambda_z
    double lmdxy;        // lambda_xy
    double alpha;        // alpha
    
    // Constructor
    GrapheneTbParams(string prefix = ""):HamParams(prefix){
        mTitle = "Graphene tight binding parameters";
    }
        
    string toString() const { 
        stringstream ss;
        ss << HamParams::toString() << ": " << endl;
        ss << mPrefix << " dtol  = " << mdtol << endl;
        ss << mPrefix << " ec    = " << ec << endl;
        ss << mPrefix << " di0   = " << di0 << endl;
        ss << mPrefix << " ti0   = " << ti0 << endl;
        ss << mPrefix << " do0   = " << do0 << endl;
        ss << mPrefix << " to0   = " << to0 << endl;
        ss << mPrefix << " doX   = " << doX << endl;
        ss << mPrefix << " lmdz  = " << lmdz << endl;
        ss << mPrefix << " lmdxy = " << lmdxy << endl; 
        ss << mPrefix << " alpha = " << alpha << endl;

        return ss.str(); 
    };   

};
/** 
 * Tight binding Hamiltonian and m for graphene.
 */
class GrapheneTbHam: public cxham{
// Methods
public:
    GrapheneTbHam(const GrapheneTbParams& p){
        mhp = shared_ptr<HamParams> (new GrapheneTbParams(p));
    };
    virtual ~GrapheneTbHam(){};

    //!< Generate Hamiltonian between two atoms.
    virtual cxmat genTwoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj);    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat genTwoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj);    
};

}
}
#endif	/* GRAPHENETB_HPP */

