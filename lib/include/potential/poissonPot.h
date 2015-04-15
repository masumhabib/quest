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
    mat mat2Phi; //!< Current Potential for the whole structure
    double nT = 300; //!< Temperature
    
private:
    double nDELX; //!< minimum difference between points in X direction
    double nDELY; //!< minimum difference between points in Y direction
    double nLx; //!< Length in X direction
    double nLy; //!< Length in Y direction
    uint nNX; //!< No of points in X direction
    uint nNY; //!< No of points in Y direction
    vec vecX; //!< Grid points in X direction
    vec vecY; //!< Grid points in Y direction
    mat mat2epsilon; //!< Dielectric constants of materials in the whole structure
    mat mat2PotentialDirichlet; //!< Dirichlet Conditions in the whole structure
    mat mat2Rho; //!< Charge density for the whole structure
    spmat mat2SparseGrad2Lambda; //!< Sparse matrix d2Lambda/dx2
    mat mat2ni; //!< intrinsic career concentration for the whole structure
    mat mat2Doping; //!< Doping concentration for the whole structure
    vec vecRHS; //!< Total Right Hand Side = Rho + d2Phi0/dx2
    vec vecd2V_by_dx2; //!< Partial Right hand side - x component
    vec vecd2V_by_dy2; //!< Partial Right hand side - y component
    vec vecRho; //!< Charge density for the whole structure in vector format
    vec vecdRho_dV; //!< Changable term in Left hand side for each iteration along 0 diagonal
    //vec RHS;
    
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
    
    vec getPotentialSliceAlongZ( double DistanceX );
    vec getPotentialSliceAlongX( double DepthZ );
    
    void helperDoping( double x1, double x2, double y1, double y2, double Doping, mat &mat2Doping );
    
    virtual string  toString() const;

private:
    double getRho( double xi, double yj, double Potential );
    void calculateh( double xi, double yj, double &hxMinus, double &hxPlus, double &hyMinus, double &hyPlus );
    
    
};
    
}
}

#endif	/* POISSONPOT_H */