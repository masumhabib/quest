/* 
 * File:   genKpAtoms.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 27, 2014, 12:33 AM
 * 
 * k.p grid points generator.
 * 
 */

#include "genKpAtoms.h"

Atoms genKpAtoms(double Lx, double Ly, double ax, double ay, 
        const vector<Atom>& periodicTable){
    MatGrid xy(-Lx/2, Lx/2, ax, -Ly/2, Ly/2, ay);
    
    // calculate x, y and z coordinates of the atoms
    int na = xy.Nx()*xy.Ny();   // total number of atoms
    mat xyz(na, 3);                 // xyz coordinate of atoms
    xyz.col(spacevec::X) = xy.X();
    xyz.col(spacevec::Y) = xy.Y();
    xyz.col(spacevec::Z).zeros();
    
    // prepare atomId list containing atomic number of a fake atom 'D'.
    icolvec atomId(na);
    atomId.fill(periodicTable[0].ia);
    
    // lattice vector
    lvec lv;
    lv.a1(spacevec::X) = xy.maxx()-xy.minx()+ax;
    lv.a2(spacevec::Y) = xy.maxy()-xy.miny()+ay;
    
    Atoms kpAtoms(atomId, xyz, lv, periodicTable);
    
    return kpAtoms;
}

