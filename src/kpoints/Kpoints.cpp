/* 
 * File:   Kpoints.cpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 * 
 * Created on February 17, 2014, 12:21 AM
 */

#include "Kpoints.h"

namespace qmicad{
    
void Kpoints::addKPoint(const point& p){
    // Read the k-points                            
    row newk(2);
    newk << p.get<0>() << p.get<1>();
    // insert new point
    mk.insert_rows(mk.n_rows,newk); 
}


void Kpoints::addKLine(const point& start, const point& end, double dk){
    // Get the starting ending and num of new k-points
    double kxs = start.get<0>();
    double kys = start.get<1>();
    
    double kxe = end.get<0>();
    double kye = end.get<1>();
    
    
    uint nknew = 1;
    if (abs(dk) > 0){
        double dlk = sqrt(pow(kxs - kxe, 2) + pow(kys - kye,2));
        nknew = uint(dlk/dk) + 1;
    }
    
    // Generate k points along the line                             
    mat newk(nknew,2);
    // generate kx and ky
    newk.col(X) = linspace<col>(kxs, kxe, nknew);
    newk.col(Y) = linspace<col>(kys, kye, nknew);
    // insert new k-points
    mk.insert_rows(mk.n_rows,newk);
}

void Kpoints::addKRect(const point& lb, const point& rt, double dk){
    double kxmin = lb.get<0>();
    double kxmax = rt.get<0>();
    double kymin = lb.get<1>();
    double kymax = rt.get<1>();

    uint nkx = 1;
    uint nky = 1;
    
    if (abs(dk) > 0){
        nkx = uint(abs(kxmax - kxmin)/dk) + 1;
        nky = uint(abs(kymax - kymin)/dk) + 1;
    }

    // Generate k mesh grid    
    row newkx(nkx);
    col newky(nky);
    // generate kx and ky
    newkx = linspace<row>(kxmin, kxmax, nkx);
    newky = linspace<col>(kymin, kymax, nky);
    // insert the new k-points
    for(int ikx = 0; ikx < nkx; ++ikx){
        mat newk(nky,2);

        newk.col(X).fill(newkx(ikx));
        newk.col(Y) = newky;
        mk.insert_rows(mk.n_rows,newk);
    }

}

    
}