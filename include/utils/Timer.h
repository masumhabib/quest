/* 
 * File:   Timer.h
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on February 11, 2014, 12:43 AM
 */

#ifndef PYTIMER_H
#define	PYTIMER_H

#include "maths/arma.hpp"
#include "utils/std.hpp"
#include "utils/stringutils.h"
#include "utils/Printable.hpp"

namespace qmicad{
namespace python{

using maths::armadillo::wall_clock;
using utils::stds::string;
using utils::Printable;
using utils::strings::ttos;


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
    
protected:
    wall_clock clock;
    double interval;
};


}
}
#endif	/* PYTIMER_H */

