/* 
 * File:   PyTISurfHam.h
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on February 13, 2014, 10:37 AM
 */

#ifndef PYTISURFHAM_H
#define	PYTISURFHAM_H

#include "PyAtomicStruct.h"
#include "../qm/kp/tikp.h"

namespace qmicad{
namespace python{

struct PyTISurfKpParams :public TISurfKpParams{

    friend ostream& operator << (ostream & out, const PyTISurfKpParams &p){
        out << p.toString();
        return out;
    }

};


class PyTISurfHam: public TISurfHam{
public:
    PyTISurfHam(const PyTISurfKpParams &tikpp):TISurfHam(tikpp){
    }
    
    void setSize(uint nbi, uint nbj = 2){
        TISurfHam::setSize(nbi, nbj);
    }
    
    void genDiagBlock(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
        uint ib)
    {
        TISurfHam::genDiagBlock(bi, bj, ib);
    }
    
    void genLowDiagBlock(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
        uint ib)
    {
        TISurfHam::genLowDiagBlock(bi, bj, ib);
    }
    
    void generate(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
                          uint ib, uint jb)
    {
        TISurfHam::generate(bi, bj, ib, jb);
    }

    
    
};

}
}
#endif	/* PYTISURFHAM_H */

