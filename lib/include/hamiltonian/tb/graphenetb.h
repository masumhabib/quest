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

class GrapheneTbParams: public cxhamparams{
public:
    // Constructor
    GrapheneTbParams( string prefix = "");        
    virtual string toString() const;  
    
    void ec(double ec) { mec = ec; update(); }
    double ec() const { return mec; }
    void di0(double di0) { mdi0 = di0; update(); }
    double di0() const { return mdi0; }
    void ti0(double do0) { mti0 = do0; update(); }
    double ti0() const { return mti0; }
    void do0(double do0) { mdo0 = do0; update(); }
    double do0() const { return mdo0; }
    void to0(double to0) { mto0 = to0; update(); }
    double to0() const { return mto0; }
    void doX(double doX) { mdoX = doX; update(); }
    double doX() const { return mdoX; }
    void lmdz(double lmdz) { mlmdz = lmdz; update(); }
    double lmdz() const { return mlmdz; }
    void lmdxy(double lmdxy) { mlmdxy = lmdxy; update(); }
    double lmdxy() const { return mlmdxy; }
    void alpha(double alpha) { malpha = alpha; update(); }
    double alpha() const { return malpha; }

    //!< Generate Hamiltonian between two atoms.
    virtual cxmat twoAtomHam(const AtomicStruct& atomi, const AtomicStruct& atomj) const;    
    //!< Generate overlap matrix between two atoms.
    virtual cxmat twoAtomOvl(const AtomicStruct& atomi, const AtomicStruct& atomj) const;        
    
private:
    //!< Default parameters.
    void setDefaultParams(){
        mdtol   = 1E-3;
        mec     = 0;
        mdi0    = 1.42;
        mti0    = 3.16;
        mdo0    = 3.35;
        mto0    = 0.39;
        mlmdz   = 0.6;
        mlmdxy  = 1.7;
        malpha  = 1.65;
        mdoX    = 6.0;
        mortho = true;
        
        mpt.add(6, "C", 1, 1);
    };
    //!< Updates internal state.
    virtual void update();

private:    
    // Tight binding parameters
    double mec;            // onside energy
    double mdi0;           // in-plane C-C bond length
    double mti0;           // in-plane C-C hopping
    double mdo0;           // out-of-plane C-C distance
    double mto0;           // out-of-plane C-C hopping
    double mdoX;           // out-of-plane C-C neighbor cut-off distance in
                          // units of 1 C-C bond length
    // Inter layer all nearest neighbors
    // See PRL 109, 236604 (2012)
    double mlmdz;         // lambda_z
    double mlmdxy;        // lambda_xy
    double malpha;        // alpha
    
};

}
}
#endif	/* GRAPHENETB_HPP */

