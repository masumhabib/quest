/* 
 * File:   Timer.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 11, 2014, 12:43 AM
 */

#ifndef PYTIMER_H
#define	PYTIMER_H

#include "../maths/arma.hpp"
#include "../utils/std.hpp"
#include "../string/stringutils.h"
#include "Printable.hpp"

namespace qmicad{
namespace python{
using maths::armadillo::wall_clock;
using namespace utils::stds;
using utils::Printable;
using utils::ttos;


class Timer: public Printable{
public:
    
    Timer(const string &prefix = ""):Printable(" " + prefix){
    }
    
    void tic(){
        clock.tic();
    }
    
    double toc(){
        interval = clock.toc();
        return interval;
    }
    
    virtual string toString() const {
        return " Runtime: " + ttos(interval) + ".";
    }
    
//    friend ostream& operator << (ostream & out, const PyTimer &p){
//        out << p.toString();
//        return out;
//    }
protected:
    wall_clock clock;
    double interval;
};


}
}
#endif	/* PYTIMER_H */

