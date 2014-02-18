/* 
 * File:   Kpoints.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 17, 2014, 12:21 AM
 */

#ifndef KPOINTS_H
#define	KPOINTS_H

#include "../maths/arma.hpp"
#include "../maths/svec.h"
#include "../maths/geometry.hpp"

namespace qmicad{
using namespace maths::spvec;
using namespace maths::armadillo;
using maths::geometry::point;

class Kpoints {
public:

    void addKPoint(const point &p);
    void addKLine(const point &start, const point &end, double dk);
    void addKRect(const point &lb, const point &rt, double dk);

protected:
    mat mk;
};

}

#endif	/* KPOINTS_H */

