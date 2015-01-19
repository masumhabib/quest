/* 
 * File:   potential.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 27, 2014, 3:47 PM
 */

#ifndef POTENTIAL_H
#define	POTENTIAL_H

#include "potential/terminal.h"

#include "maths/geometry.hpp"
#include "maths/svec.h"

#include "atoms/AtomicStruct.h"

#include "utils/std.hpp"
#include "utils/stringutils.h"


namespace qmicad{
namespace potential{

using utils::strings::itos;
using namespace atoms;
using namespace utils::stds;
using namespace maths::spvec;

/**
 * Potential class handles electrostatic potential.
 */

class Potential:public Printable{

public:
protected:
    AtomicStruct::ptr   ma;     //!< Atomistic geometry of the device.
    vec                 mV;     //!< Electrostatic potential.
    vec                 mRho;   //!< Electric charge.
    vector<contact>     ms;     //!< Source contact.
    vector<contact>     md;     //!< Drain contact.
    vector<gate>        mg;     //!< Gates.  
    
public:
    //<!< Constructor.
    Potential(AtomicStruct::ptr atoms = AtomicStruct::ptr(), const string &prefix = "");
    //!< Convert atomic potential to orbital potential.
    shared_ptr<vec> toOrbPot(span s = span::all);
    vec toOrbPot(uint start, uint end);
    double Vatom(uint ia);
    void Vatom(uint ia, double V);
    //!< Calculate electrostatic potential.
    virtual void compute(){};
    //!< Convert to string for cout.
    virtual string toString() const;
    //!< Export geometry to SVG file.
    virtual void exportSvg(const string &path);
    //!< Export potential to text file.
    virtual void exportPotential(const string &path);
    //!<
    void addSource(const squadrilateral &sq);
    void addDrain(const squadrilateral &sq);
    void addGate(const squadrilateral &sq);
    
    //!< Sets the gate voltage.
    void VG(int ig, double VG);
    //!< Sets electrostatic potential at the drain.
    void VD(double VD);
    //!< Sets electrostatic potential at the source.
    void VS(double VS);
    
    uint NG() const { return mg.size(); };
    
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
}

#endif	/* POTENTIAL_H */

