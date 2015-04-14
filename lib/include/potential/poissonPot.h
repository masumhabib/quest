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
    poissonPot( vec Lx, vec Ly);
    poissonPot( double nLx, double nLy, double DELX, double DELY );
    void setMaterialEps( double x1, double x2, double y1, double y2, double epsilonR );
    void setMaterialni( double x1, double x2, double y1, double y2, double ni );
    void setDoping( double x1, double x2, double y1, double y2, double Nad );
    
    void setPotentialDirichlet( double x1, double x2, double y1, double y2, double V );
    void setPhiByFroce( double x1, double x2, double y1, double y2, double V );
    void generateGrad2LambdaMatrix( );
    mat setInitialGuess( );
    mat calculateLambdaSingleIteration( );
    double getRho( double xi, double yj, double Potential );
    vec getPotentialSliceAlongZ( double DistanceX );
    vec getPotentialSliceAlongX( double DepthZ );
    
    void helperDoping( double x1, double x2, double y1, double y2, double Doping, mat &mat2Doping );
    
    
    
    virtual string  toString() const;
private:
    void calculateh( double xi, double yj, double &hxMinus, double &hxPlus, double &hyMinus, double &hyPlus );
    
};
    
}
}

#endif	/* POISSONPOT_H */