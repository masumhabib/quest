/* 
 * File:   PyGrapheneKpHam.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 10, 2014, 12:31 AM
 */

#ifndef PYGRAPHENEKPHAM_H
#define	PYGRAPHENEKPHAM_H

#include "PyAtomicStruct.h"
#include "../qm/kp/graphenekp.h"
#include "../utils/std.hpp"

namespace qmicad{
namespace python{

struct PyGrapheneKpParams :public GrapheneKpParams{

    friend ostream& operator << (ostream & out, const PyGrapheneKpParams &p){
        out << p.toString();
        return out;
    }

};

class PyGrapheneKpHam:public GrapheneKpHam{
public:
    PyGrapheneKpHam(const PyGrapheneKpParams &gkpp):GrapheneKpHam(gkpp){
    }

    void setSize(uint nbi, uint nbj){
        GrapheneKpHam::setSize(nbi, nbj);
    }
    
    void setSizeForNegf(uint nbi){
        GrapheneKpHam::setSizeForNegf(nbi);
    }

    void setSizeForBand(uint nbi){
        GrapheneKpHam::setSizeForBand(nbi);
    }
    
    void genDiagBlock(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
        uint ib)
    {
        GrapheneKpHam::genDiagBlock(bi, bj, ib);
    }
    
    void genLowDiagBlock(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
        uint ib)
    {
        GrapheneKpHam::genLowDiagBlock(bi, bj, ib);
    }

    void genNearestNeigh(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
        uint ib)
    {
        GrapheneKpHam::genNearestNeigh(bi, bj, ib);
    }
    
    void generate(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
                          uint ib, uint jb)
    {
        GrapheneKpHam::generate(bi, bj, ib, jb);
    }
    
};

}
}
#endif	/* PYGRAPHENEKPHAM_H */

