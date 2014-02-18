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
        Printable(" " + prefix), ma(atoms), ms(source), md(drain), mg(gates)
{
    mV.set_size(ma.NumOfAtoms());
    mV.zeros();
    ms.Title("Source");
    md.Title("Drain");
    
    ms.Prefix(ms.Prefix() + mPrefix);
    md.Prefix(md.Prefix() + mPrefix);
    
    for (int it = 0; it < mg.size(); ++it){
        mg[it].Title("Gate # " + itos(it));
        mg[it].Prefix(mg[it].Prefix() + mPrefix);
    }    
}

Potential::Potential(const AtomicStruct &atoms, const string &prefix):
        Printable(" " + prefix), ma(atoms)
{
    mV.set_size(ma.NumOfAtoms());
    mV.zeros();
}


string Potential::toString() const{
    stringstream ss;
    ss << Printable::toString() << ":" << endl;
    ss << mPrefix << ms << endl;
    for (vector<gate>::const_iterator it = mg.begin(); it != mg.end(); ++it){
        ss << mPrefix << *it << endl;
    }
    ss << mPrefix << md << endl;
   
    return ss.str();
}

// Convert atomic potential to orbital potential
shared_ptr<vec> Potential::toOrbPot(span s){
    AtomicStruct a = ma(s);
    int no = a.NumOfOrbitals();
    int na = a.NumOfAtoms();
    shared_ptr<vec> pV = make_shared<vec>(no, fill::zeros);
    vec &V = *pV;

    int io = 0;
    // loop through the atoms we are interested in
    for (int ia = 0; ia < na; ++ia){
        int n = a.AtomAt(ia).no;
        // fill the orbitals of atom ia with its potential
        V(span(io,io+n-1)).fill(mV(ia + s.a));
        io += n;
    }  

    return pV;
}

shared_ptr<vec> Potential::toOrbPot(uint start, uint end){
    return toOrbPot(span(start, end));
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
    Va.col(coord::X) = ma.X();
    Va.col(coord::Y) = ma.Y();
    Va.col(coord::Z) = ma.Z();
    Va.col(coord::Z+1) = mV;
    
    ofstream potf(path.c_str());
    potf << Va << endl;
}

void Potential::addSource(const squadrilateral &sq){
    // source
    ms = contact(sq.lb, sq.rb, sq.rt, sq.lt);
    ms.Title("Source");
    ms.Prefix(ms.Prefix() + mPrefix);
}

void Potential::addDrain(const squadrilateral &sq) {
    // drain
    md = contact(sq.lb, sq.rb, sq.rt, sq.lt);
    md.Title("Drain");
    md.Prefix(md.Prefix() + mPrefix);
}

void Potential::addGate(const squadrilateral& sq){
    mg.push_back(gate(sq.lb, sq.rb, sq.rt, sq.lt)); 
    int it = mg.size() - 1;
    mg[it].Title("Gate # " + itos(it+1));
    mg[it].Prefix(mg[it].Prefix() + mPrefix);
}

void Potential::VG(int ig, double VG){
    mg[ig].V = VG;
}

void Potential::VD(double VD){
    md.V = VD;
}
void Potential::VS(double VS){
    ms.V = VS;
}

}

