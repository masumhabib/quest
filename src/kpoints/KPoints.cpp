/* 
 * File:   KPoints.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on February 17, 2014, 12:21 AM
 */

#include "kpoints/KPoints.h"

namespace qmicad{
namespace kpoints{

KPoints::KPoints(const string& prefix):Printable(" " + prefix){
    
}
    
void KPoints::addKPoint(const point& p){
    // Read the k-points                            
    row newk(3);
    newk << p.get<0>() << p.get<1>() << 0.0;
    // insert new point
    mk.insert_rows(mk.n_rows,newk); 
}


void KPoints::addKLine(const point& start, const point& end, double dk){
    double dkx = start.get<0>() - end.get<0>();
    double dky = start.get<1>() - end.get<1>();
    
    // Compute num of new k-points    
    uint nk = 1;
    if (abs(dk) > 0){
        double dlk = sqrt(pow(dkx, 2) + pow(dky,2));
        nk = uint(dlk/dk) + 1;
    }
    
    addKLine(start, end, nk);
}

void KPoints::addKLine(const point& start, const point& end, uint nk){
    // Get the starting ending and num of new k-points
    double kxs = start.get<0>();
    double kys = start.get<1>();
    
    double kxe = end.get<0>();
    double kye = end.get<1>();
        
    // Generate k points along the line                             
    mat newk(nk,3);
    // generate kx and ky
    newk.col(coord::X) = linspace<col>(kxs, kxe, nk);
    newk.col(coord::Y) = linspace<col>(kys, kye, nk);
    newk.col(coord::Z).zeros();
    
    // insert new k-points
    mk.insert_rows(mk.n_rows,newk);
}

void KPoints::addKRect(const point& lb, const point& rt, double dkx, double dky){
    double delkx = rt.get<0>() - lb.get<0>();
    double delky = rt.get<1>() - lb.get<1>();

    uint nkx = 1;
    uint nky = 1;
    
    if (abs(dkx) > 0){
        nkx = uint(abs(delkx/dkx)) + 1;
    }

    if (abs(dky) > 0){
        nky = uint(abs(delky/dky)) + 1;
    }

    addKRect(lb, rt, nkx, nky);
}

void KPoints::addKRect(const point& lb, const point& rt, uint nkx, uint nky){
    double kxmin = lb.get<0>();
    double kxmax = rt.get<0>();
    double kymin = lb.get<1>();
    double kymax = rt.get<1>();

    // Generate k mesh grid    
    row newkx(nkx);
    col newky(nky);
    // generate kx and ky
    newkx = linspace<row>(kxmin, kxmax, nkx);
    newky = linspace<col>(kymin, kymax, nky);
    // insert the new k-points
    for(uint ikx = 0; ikx < nkx; ++ikx){
        mat newk(nky,3);

        newk.col(coord::X).fill(newkx(ikx));
        newk.col(coord::Y) = newky;
        newk.col(coord::Z).zeros();         
        mk.insert_rows(mk.n_rows,newk);
    }
}

shared_ptr<mat> KPoints::kp(){
    shared_ptr<mat> k(&mk, NullDeleter());
    return k;
}
    
}
}




