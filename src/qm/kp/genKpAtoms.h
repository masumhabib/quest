/* 
 * File:   kpatoms.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 27, 2014, 12:32 AM
 * 
 * k.p grid points generator.
 * 
 */

#ifndef KPATOMS_H
#define	KPATOMS_H

#include "../../grid/grid.hpp"
#include "../../atoms/Atoms.h"
#include "../../atoms/Lattice.h"
#include "../../utils/svec.h"

/*
 * Generates discretized k.p lattice points and associates them to 
 * fake atoms.
 */
Atoms genKpAtoms(double Lx, double Ly, double ax, double ay, 
                 const vector<Atom>& periodicTable);

#endif	/* KPATOMS_H */

