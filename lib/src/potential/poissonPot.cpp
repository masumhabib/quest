/*
 * File:   poissonPot.cpp
 * Author: Mirza Elahi <mirza.monzur@gmail.com>
 *
 * Created on March 29, 2015, 6:34 AM
 */

#include "potential/poissonPot.h"

namespace qmicad{
namespace potential{
    
    
poissonPot::poissonPot(vec vecLx, vec vecLy) {
    using namespace utils::stds;
    this->mTitle = "Poisson Voltage Profile";
    
    this->vecX = vecLx;
    this->vecY = vecLy;
    this->nNX = vecLx.n_elem;
    this->nNY = vecLy.n_elem;
    
    this->nLx = vecLx( this->nNX - 1 );
    this->nLy = vecLy( this->nNY - 1 );
    
    
    
}
    
poissonPot::poissonPot(double nLx, double nLy, double DELX, double DELY) {
    
}
    
void poissonPot::setMaterialEps( double x1, double x2, double y1, double y2, double epsilonR ){
    
    
}
    
void poissonPot::setMaterialni( double x1, double x2, double y1, double y2, double ni ){
    
}
    
void poissonPot::setDoping( double x1, double x2, double y1, double y2, double Nad ){
    
}

void poissonPot::setPotentialDirichlet( double x1, double x2, double y1, double y2, double V ){
    
}

void poissonPot::setPhiByFroce( double x1, double x2, double y1, double y2, double V ){
    
}

void poissonPot::generateGrad2LambdaMatrix( ){
    
}

mat poissonPot::setInitialGuess( ){
    mat A;
    return A;
}

mat poissonPot::calculateLambdaSingleIteration( ){
    mat A;
    return A;
}

double poissonPot::getRho( double xi, double yj, double Potential ){
    
    return 5.0;
}

vec poissonPot::getPotentialSliceAlongZ( double DistanceX ){
    vec A;
    return A;
}

vec poissonPot::getPotentialSliceAlongX( double DepthZ ){
    vec A;
    return A;
}

void poissonPot::helperDoping( double x1, double x2, double y1, double y2, double Doping, mat &mat2Doping ){
    
}

void calculateh( double xi, double yj, double &hxMinus, double &hxPlus, double &hyMinus, double &hyPlus ){
    
}
    

string poissonPot::toString() const
{
    stringstream ss;
    ss << Printable::toString() << ":" << endl;
    ss << mPrefix << " Number of gates: " << NG() << endl;
//    ss << mPrefix << " Number of linear regions: " << NLR() << endl;
//    for (vector<contact>::const_iterator it = ms.begin(); it != ms.end(); ++it){
//        ss << mPrefix << *it << endl;
//    }
//    
//    for (vector<gate>::const_iterator it = mg.begin(); it != mg.end(); ++it){
//        ss << mPrefix << *it << endl;
//    }
//    for (vector<linear_region>::const_iterator it = mlr.begin(); it != mlr.end(); ++it){
//        ss << mPrefix << *it << endl;
//    }
//    
//    for (vector<contact>::const_iterator it = md.begin(); it != md.end(); ++it){
//        ss << mPrefix << *it << endl;
//    }
    
    return ss.str();
}
    

    
    
    
    
}
}