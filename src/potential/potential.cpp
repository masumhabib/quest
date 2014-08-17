/* 
 * File:   potential.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 27, 2014, 3:47 PM
 */

#include "potential/potential.h"

namespace qmicad{
namespace potential{

using namespace std;

Potential::Potential(AtomicStruct::ptr atoms, const string &prefix):
        Printable(" " + prefix), ma(atoms)
{
    if (ma){
        mV.set_size(ma->NumOfAtoms());
        mV.zeros();

        mRho.set_size(ma->NumOfAtoms());
        mRho.zeros();
    }

}


string Potential::toString() const{
    stringstream ss;
    ss << Printable::toString() << ":" << endl;
    for (vector<contact>::const_iterator it = ms.begin(); it != ms.end(); ++it){
        ss << mPrefix << *it << endl;
    }

    for (vector<gate>::const_iterator it = mg.begin(); it != mg.end(); ++it){
        ss << mPrefix << *it << endl;
    }

    for (vector<contact>::const_iterator it = md.begin(); it != md.end(); ++it){
        ss << mPrefix << *it;
    }
   
    return ss.str();
}

// Convert atomic potential to orbital potential
shared_ptr<vec> Potential::toOrbPot(span s){
    AtomicStruct a = (*ma)(s);
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

vec Potential::toOrbPot(uint start, uint end){
    return *toOrbPot(span(start, end));
}

double Potential::Vatom(uint ia){
    return mV(ia);
}

void Potential::Vatom(uint ia, double V){
    mV(ia) = V;
}

void Potential::exportSvg(const string& path){
    using namespace std;
    using namespace bg;
    ofstream svg(path.c_str());
    
    svg_mapper<point> mapper(svg, 1024, 768);
    
    // add polygons to mapper
    for (int it = 0; it < ms.size(); ++it){
        mapper.add(ms[it].geom);
    }
    
    for (int it = 0; it < mg.size(); ++it){
        mapper.add(mg[it].geom);
    }
    
    for (int it = 0; it < md.size(); ++it){
        mapper.add(md[it].geom);
    }
    
    // draw polygons on mapper
    for (int it = 0; it < ms.size(); ++it){
        mapper.map(ms[it].geom, "fill-opacity:0.4;fill:rgb(10,10,255);stroke:rgb(10,10,255);stroke-width:2");
    }
    
    for (int it = 0; it < mg.size(); ++it){
        mapper.map(mg[it].geom, "fill-opacity:0.4;fill:rgb(204,10,204);stroke:rgb(204,10,204);stroke-width:2");
    }
    
    for (int it = 0; it < md.size(); ++it){
        mapper.map(md[it].geom, "fill-opacity:0.4;fill:rgb(255,10,10);stroke:rgb(255,10,10);stroke-width:2");
    }

}

void Potential::exportPotential(const string& path){
    int na = ma->NumOfAtoms();
    mat Va(na, 4);
    Va.cols(coord::X, coord::Z) = ma->XYZ();
    Va.col(coord::Z+1) = mV;
    
    ofstream potf(path.c_str());
    potf << Va << endl;
}

void Potential::addSource(const squadrilateral &sq){
    // source
    if (ms.empty()){
        ms.push_back(contact(sq.lb, sq.rb, sq.rt, sq.lt));
    }else{
        ms[0] = contact(sq.lb, sq.rb, sq.rt, sq.lt);
    }
    ms[0].Title("Source");
    ms[0].Prefix(ms[0].Prefix() + mPrefix);
}

void Potential::addDrain(const squadrilateral &sq) {
    // drain
    if(md.empty()){
        md.push_back(contact(sq.lb, sq.rb, sq.rt, sq.lt));
    }else{
        md[0] = contact(sq.lb, sq.rb, sq.rt, sq.lt);
    }
    md[0].Title("Drain");
    md[0].Prefix(md[0].Prefix() + mPrefix);
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
    md[0].V = VD;
}
void Potential::VS(double VS){
    ms[0].V = VS;
}

}
}


    

