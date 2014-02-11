/* 
 * File:   potential.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 27, 2014, 3:47 PM
 */

#include "potential.h"

namespace qmicad{
using namespace std;
// Constructor
Potential::Potential(const AtomicStruct &atoms, const contact &source, 
const contact &drain, const vector<gate> &gates, const string &prefix):
Printable("" + prefix), ma(atoms), ms(source), md(drain), mg(gates){
    mV.set_size(ma.NumOfAtoms());
    mV.zeros();
    ms.Title("Source");
    md.Title("Drain");
    
    ms.Prefix(" " + prefix);
    md.Prefix(" " + prefix);
    
    for (int it = 0; it < mg.size(); ++it){
        mg[it].Title("Gate # " + itos(it));
        mg[it].Prefix(" " + prefix);
    }    
}


string Potential::toString() const{
    stringstream ss;
    ss << mPrefix << ms << endl;
    for (vector<gate>::const_iterator it = mg.begin(); it != mg.end(); ++it){
        ss << mPrefix << *it << endl;
    }
    ss << mPrefix << md << endl;
   
    return ss.str();
}

// Convert atomic potential to orbital potential
vec Potential::toOrbPot(span s){
    AtomicStruct a = ma(s);
    int no = a.NumOfOrbitals();
    int na = a.NumOfAtoms();
    vec V(no, fill::zeros);

    int io = 0;
    // loop through the atoms we are interested in
    for (int ia = 0; ia < na; ++ia){
        int n = a.AtomAt(ia).no;
        // fill the orbitals of atom ia with its potential
        V(span(io,io+n-1)).fill(mV(ia + s.a));
        io += n;
    }  

    return V;
}

void Potential::exportSvg(const string& path){
    using namespace std;
    using namespace bg;
    ofstream svg(path.c_str());
    
    svg_mapper<point> mapper(svg, 1024, 768);
    
    // add polygons to mapper
    mapper.add(ms.geom);
    for (int it = 0; it < mg.size(); ++it){
        mapper.add(mg[it].geom);
    }
    mapper.add(md.geom);
    
    // draw polygons on mapper
    mapper.map(ms.geom, "fill-opacity:0.4;fill:rgb(10,10,255);stroke:rgb(10,10,255);stroke-width:2");
    for (int it = 0; it < mg.size(); ++it){
        mapper.map(mg[it].geom, "fill-opacity:0.4;fill:rgb(204,10,204);stroke:rgb(204,10,204);stroke-width:2");
    }
    mapper.map(md.geom, "fill-opacity:0.4;fill:rgb(255,10,10);stroke:rgb(255,10,10);stroke-width:2");

}

void Potential::exportPotential(const string& path){
    int na = ma.NumOfAtoms();
    mat Va(na, 4);
    Va.col(spacevec::X) = ma.X();
    Va.col(spacevec::Y) = ma.Y();
    Va.col(spacevec::Z) = ma.Z();
    Va.col(spacevec::Z+1) = mV;
    
    ofstream potf(path.c_str());
    potf << Va << endl;
}

}

