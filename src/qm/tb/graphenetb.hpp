/* 
 * File:   graphenetb.hpp
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

#include "../../atoms/Atoms.h"
#include "../hamParams.hpp"

namespace qmicad{
using arma::mat;
    
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
    
    // our periodic table
    vector<Atom> PeriodicTable;
    
    
    // Constructor
    GrapheneTbParams(string prefix = "") :HamParams(prefix){
        mTitle = "Graphene";
        // default parameters
        
        PeriodicTable.push_back(Atom(6,  "C",  1, 1));

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
/* Tight binding logic for graphene*/
class GrapheneTbHamGen: public HamGenerator<mat>{
// Fields
protected:    
    GrapheneTbParams p;
// Methods
public:
    GrapheneTbHamGen(const GrapheneTbParams& p){this->p = p;};
    virtual ~GrapheneTbHamGen(){};
 
    virtual mat operator()(const Atoms& atomi, const Atoms& atomj){
        float d, dx, dy, dz;
        int noi = atomi.NumOfOrbitals();
        int noj = atomj.NumOfOrbitals();
        
        mat hmat =  zeros<mat>(noi, noj);
        
        // calculate distance between atom i and atom j
        dx = atomi.X(0) - atomj.X(0);
        dy = atomi.Y(0) - atomj.Y(0);           
        dz = atomi.Z(0) - atomj.Z(0);
        d = sqrt(dx*dx + dy*dy + dz*dz);

        // Assign the the matrix elements based on the distance between 
        // the atoms
        // C-C
        if (atomi.Symbol(0) == "C" && atomj.Symbol(0) == "C"){

            dz = abs(dz); // inter-plane distance

            // site energy
            if (d <= p.dtol){ 
                hmat(0,0) = p.ec;
            // in-plane first nearest neighbor
            }else if (dz <= p.dtol && abs(d - p.di0) <= p.dtol){
                hmat(0,0) = -p.ti0;
            /*// out-of-plane first nearest neighbors
            }else if (abs(delz - do0cc) <= dtol && abs(d - do0cc) <= dtol){
                hmat(ia,ja) = -to0cc;
            }*/
            // out-of-plane some nearest neighbors
            }else if (abs(dz - p.do0) <= p.dtol && abs(d - p.do0) <= (p.doX*p.di0 + p.dtol)){
                // PRB 84, 195421 (2011)
                // DBG hmat(ia,ja) = -to0cc*exp(-3.0*(d - do0cc));

                // PRL 109, 236604 (2012)
                double dxy = sqrt(dx*dx + dy*dy);
                hmat(0,0) = -p.to0*exp(-(d - p.do0)/p.lmdz)*exp(-pow(dxy/p.lmdxy, p.alpha));
            }
        }

    };
    
private:

};

}
#endif	/* GRAPHENETB_HPP */

