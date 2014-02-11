/* 
 * File:   potential.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 27, 2014, 3:47 PM
 */

#ifndef POTENTIAL_H
#define	POTENTIAL_H

#include "terminal.h"

#include "../string/stringutils.h"
#include "../maths/geometry.hpp"
#include "../atoms/AtomicStruct.h"

namespace qmicad{
using utils::itos;

/**
 * Potential class handles electrostatic potential.
 */

class Potential:public Printable{
public:
protected:
    const AtomicStruct  &ma;    //<! Atomistic geometry of the device
    contact             ms;     //<! Source contact
    contact             md;     //<! Drain contact
    vector<gate>        mg;     //<! Gates    
    vec                 mV;     //<! Electrostatic potential
    
public:
    //<! Constructor
    Potential(const AtomicStruct &atoms, const contact &source, 
            const contact &drain, const vector<gate> &gates, 
            const string &prefix = "");
    //<! Constructor
    Potential(const AtomicStruct &atoms, const string &prefix = "");
    //<! Convert atomic potential to orbital potential
    vec toOrbPot(span s = span::all);
    //<! Calculate electrostatic potential
    virtual void compute(){};
    //<! Convert to string for cout
    virtual string toString() const;
    //<! Export geometry to SVG file
    virtual void exportSvg(const string &path);
    //<! Export potential to text file
    virtual void exportPotential(const string &path);
    //<!
    void addSource(const squadrilateral &sq);
    void addDrain(const squadrilateral &sq);
    void addGate(const squadrilateral& sq);
    
protected:
    struct Contains{
        point p;
        Contains(double x, double y):p(x,y){};
        bool operator ()(Terminal &T){
            return bg::within(p, T.geom);
        }
    };
    
};

}

#endif	/* POTENTIAL_H */

