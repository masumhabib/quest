/* 
 * File:        linearPot.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
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

#include "potential/terminal.h"
#include "potential/potential.h"
#include "utils/vout.h"

#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>



namespace qmicad{
namespace potential{

using namespace utils::stds;
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
    LinearPot(AtomicStruct::ptr atoms = AtomicStruct::ptr(), const string &prefix = "");
    
    virtual string  toString() const;
    virtual void    exportSvg(const string &path);
    virtual void    compute();  
    
    void addLinearRegion(const squadrilateral& sq);
    void VLR(int ilr, double Vl, double Vr, double Vt = 0, double Vb = 0);
    
    uint NLR() const { return mlr.size(); };

};

}
}
#endif	/* LINEARPOT_H */

