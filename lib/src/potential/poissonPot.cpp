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
    
    this->nDELX = min(   vecLx.rows( 1 , this->nNX - 1 ) - vecLx.rows( 0 , this->nNX - 2 )   );
    this->nDELY = min(   vecLy.rows( 1 , this->nNY - 1 ) - vecLy.rows( 0 , this->nNY - 2 )   );
    
    this->initAuxMatrices();  
}

poissonPot::poissonPot(double nLx, double nLy, double DELX, double DELY) {
    this->vecX = linspace<double> (0, nLx, DELX);
    this->vecY = linspace<double> (0, nLy, DELY);
    this->nLx = this->vecX.max();
    this->nLy = this->vecY.max();
    this->nDELX = DELX;
    this->nDELY = DELY;
    
    this->initAuxMatrices();
}
 
void poissonPot::initAuxMatrices() {
    this->mat2epsilon            = zeros<mat>( this->nNX, this->nNY );
    this->mat2Doping             = zeros<mat>( this->nNX, this->nNY );
    this->mat2Rho                = zeros<mat>( this->nNX, this->nNY );
    this->mat2PotentialDirichlet = this->Inf * ones<mat>( this->nNX, this->nNY ); 
    this->mat2ni                 = zeros<mat>( this->nNX, this->nNY );
    this->mat2SparseGrad2Lambda  = spmat( this->nNX * this->nNX, this->nNY * this->nNY  );
    this->mat2Phi                = zeros<mat>( this->nNX, this->nNY );
    this->vecRHS                 = zeros<vec>( this->nNX * this->nNY );
    this->vecd2V_by_dx2          = zeros<vec>( this->nNX * this->nNY );
    this->vecd2V_by_dy2          = zeros<vec>( this->nNX * this->nNY );
    this->vecRho                 = zeros<vec>( this->nNX * this->nNY );
}


    
void poissonPot::setMaterialEps( double x1, double x2, double y1, double y2, double epsilonR ){
    
    uvec tempXup = find ( this->vecX > x1 );
    uvec tempXLow = find ( this->vecX <= y1 );
    uvec rowIndices;
    maths::vintersection( tempXup, tempXLow, rowIndices );
    
    
    uvec tempYup = find ( this->vecY > y1 );
    uvec tempYLow = find ( this->vecY <= y2 );
    uvec colIndices;
    maths::vintersection( tempYup, tempYLow, colIndices );
    
    this->mat2epsilon( rowIndices, colIndices ).fill( epsilonR );   
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
    ss << mPrefix << " Total Length in x direction " << this->nLx << endl;
    ss << mPrefix << " Total Length in y direction " << this->nLy << endl;
    ss << mPrefix << " Minimum difference in x direction " << this->nDELX << endl;
    ss << mPrefix << " Minimum difference in y direction " << this->nDELY << endl;
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