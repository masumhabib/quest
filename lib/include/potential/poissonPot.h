/*
 * File:   poissonPot.h
 * Author: Mirza Elahi <mirza.monzur@gmail.com>
 *
 * Created on March 29, 2015, 6:34 AM
 */

#ifndef POISSONPOT_H
#define	POISSONPOT_H

#include "potential/terminal.h"
#include "potential/potential.h"
#include "utils/vout.h"

#include <vector>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>

using namespace utils::stds;

namespace qmicad{
    namespace potential{
        
        
        class poissonPot : public Potential {
        public:
            poissonPot(double Lx, double Ly);
            
        private:
            
        };
        
    }
}

#endif	/* POISSONPOT_H */