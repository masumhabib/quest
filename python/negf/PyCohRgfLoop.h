/* 
 * File:   PyCohRgfLoop.h
 * Author: K M Masum Habib <masum.habib@virginia.edu>
 *
 * Created on August 17, 2014, 10:17 AM
 */

#ifndef PYCOHRGFLOOP_H
#define	PYCOHRGFLOOP_H

#include "negf/CohRgfLoop.h"
#include "negf/RgfResult.h"
#include "boostpython.hpp"
#include "npyarma/npyarma.h"

namespace qmicad{
namespace python{

using namespace negf;
namespace bp = boost::python;
using qmicad::python::npy2mat;
using qmicad::python::npy2col;


class PyCohRgfLoop: public CohRgfLoop{
public:
    PyCohRgfLoop(const Workers &workers, uint nb = 5, double kT = 0.0259, 
        dcmplx ieta = dcmplx(0,1E-3), bool orthogonal = true, uint nTransNeigh = 0,
        string newprefix = ""); 
    
    // For python binding
//    void            H0(bp::object H0, int ib, int ineigh);
//    void            S0(bp::object S0, int ib, int ineigh);
//    void            Hl(bp::object Hl, int ib, int ineigh);
//    void            Sl(bp::object Sl, int ib, int ineigh);    
//    void            H0(bp::object H0, int ib);
//    void            S0(bp::object S0, int ib);
//    void            Hl(bp::object Hl, int ib);
//    void            Sl(bp::object Sl, int ib);    
    void            H0(const cxmat& H0, int ib, int ineigh);
    void            S0(const cxmat& S0, int ib, int ineigh);
    void            Hl(const cxmat& Hl, int ib, int ineigh);
    void            Sl(const cxmat& Sl, int ib, int ineigh);    
    void            H0(const cxmat& H0, int ib);
    void            S0(const cxmat& S0, int ib);
    void            Hl(const cxmat& Hl, int ib);
    void            Sl(const cxmat& Sl, int ib);    

    void            V(const col& Sl, int ib);
    void            pv0(const col& pv0, int ib, int ineigh);
    void            pvl(const col& pv0, int ib, int ineigh);

};

}}

#endif	/* PYCOHRGFLOOP_H */

