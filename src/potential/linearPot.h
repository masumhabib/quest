/* 
 * File:        linearPot.h
 * Author:      K M Masum Habib<masum.habib@virginia.edu>
 * 
 * Description: 
 * Linear potential class generates linear potential profile
 * throughout the device using dirichlet boundary condition at gates and
 * Neumann boundary condition at source and drain. The gates are inherently
 * quadrilateral.
 * 
 * @FIXME Only works for rectangular gates.
 * 
 * Created on January 27, 2014, 5:19 PM
 * Last modified: 
 */

#ifndef LINEARPOT_H
#define	LINEARPOT_H

#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "terminal.h"
#include "potential.h"
#include "../string/stringutils.h"


namespace qmicad{
using utils::dtos;
/*
 * Linear terminal emulates linear voltage drop between two 
 * gates.
 */
struct LinearRegion4:public Terminal4{
    double Vl;          // left gate voltage.
    double Vr;          // right gate voltage.
    double Vt;          // top gate voltage.
    double Vb;          // bottom gate voltage.
    
    LinearRegion4(point lb, point rb, point rt, point lt, 
            double Vl = 0, double Vr = 0, double Vt = 0, double Vb = 0,
            const string &prefix = "");
    
    virtual string toString() const;
};

typedef LinearRegion4 linear_region;

/*
 * Linear potential drop model
 */

class LinearPot:public Potential{
protected:
    vector<linear_region> mlr; // linear voltage region
    
public:
    LinearPot(const AtomicStruct &atoms, const contact &source, 
            const contact &drain, const vector<gate> &gates, 
            const vector <linear_region> &linear, const string &prefix = "");
 
    LinearPot(const AtomicStruct &atoms, const string &prefix = "");
    
    virtual string  toString() const;
    virtual void    exportSvg(const string &path);
    virtual void    compute();  
    
    void addLinearRegion(const squadrilateral& sq);
    void VLR(int ilr, double Vl, double Vr, double Vt = 0, double Vb = 0);
};

}
#endif	/* LINEARPOT_H */

