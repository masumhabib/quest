/* 
 * File:   KPoints.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 17, 2014, 12:21 AM
 */

#ifndef KPOINTS_H
#define	KPOINTS_H

#include "maths/arma.hpp"
#include "maths/svec.h"
#include "maths/geometry.hpp"
#include "utils/NullDeleter.hpp"
#include "utils/vout.h"
#include "utils/Printable.hpp"

#include <boost/shared_ptr.hpp>

namespace qmicad{
namespace kpoints{

using namespace utils::stds;
using utils::NullDeleter;
using utils::Printable;

using namespace maths::spvec;
using namespace maths::armadillo;
using maths::geometry::point;

using boost::shared_ptr;

class KPoints: public Printable {
public:

    KPoints(const string& prefix = "");
    
    void            addKPoint(const point &p);
    void            addKLine(const point& start, const point& end, double dk);
    void            addKLine(const point& start, const point& end, uint nk);
    void            addKRect(const point& lb, const point& rt, double dkx, double dky);
    void            addKRect(const point& lb, const point& rt, uint nkx, uint nky);
    
    shared_ptr<mat> kp();
    uint            N() {return mk.n_rows; }
    
    virtual string toString(){
        stringstream out;
        out << endl << mk;
        return out.str();
    }

protected:
    mat mk;
};

}
}
#endif	/* KPOINTS_H */

