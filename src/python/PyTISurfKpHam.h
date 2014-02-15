/* 
 * File:   PyTISurfKpHam.h
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


class PyTISurfKpHam: public TISurfKpHam{
public:
    PyTISurfKpHam(const PyTISurfKpParams &tikpp):TISurfKpHam(tikpp){
    }
    
    void setSize(uint nbi, uint nbj = 2){
        TISurfKpHam::setSize(nbi, nbj);
    }
    
    void genDiagBlock(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
        uint ib)
    {
        TISurfKpHam::genDiagBlock(bi, bj, ib);
    }
    
    void genLowDiagBlock(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
        uint ib)
    {
        TISurfKpHam::genLowDiagBlock(bi, bj, ib);
    }
    
    void generate(const PyAtomicStruct &bi, const PyAtomicStruct &bj, 
                          uint ib, uint jb)
    {
        TISurfKpHam::generate(bi, bj, ib, jb);
    }

    
    
};

}
}
#endif	/* PYTISURFHAM_H */

