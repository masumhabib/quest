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
    this->vecdRho_dV             = zeros<vec>( this->nNX * this->nNY );
}
  
void poissonPot::setMaterialEps( double x1, double x2, double y1, double y2, double epsilonR ){
    
    uwcol rowIndices;
    uwcol colIndices;
    this->helperGenRowColIndices( x1, x2, y1, y2, rowIndices, colIndices);
    this->mat2epsilon( rowIndices, colIndices ).fill( epsilonR );   
}

void poissonPot::setMaterialni( double x1, double x2, double y1, double y2, double ni ){
    uwcol rowIndices;
    uwcol colIndices;
    this->helperGenRowColIndices( x1, x2, y1, y2, rowIndices, colIndices);
    this->mat2ni( rowIndices, colIndices ).fill( ni ); 
}
    
void poissonPot::setDoping( double x1, double x2, double y1, double y2, double Nad ){
    uwcol rowIndices;
    uwcol colIndices;
    this->helperGenRowColIndices( x1, x2, y1, y2, rowIndices, colIndices);
    this->mat2Doping( rowIndices, colIndices ).fill( Nad ); 
}

void poissonPot::setPotentialDirichlet( double x1, double x2, double y1, double y2, double V ){
    uwcol rowIndices;
    uwcol colIndices;
    this->helperGenRowColIndices( x1, x2, y1, y2, rowIndices, colIndices);
    this->mat2PotentialDirichlet( rowIndices, colIndices ).fill( V ); 
}

void poissonPot::setPhiByFroce( double x1, double x2, double y1, double y2, double V ){
    uwcol rowIndices;
    uwcol colIndices;
    this->helperGenRowColIndices( x1, x2, y1, y2, rowIndices, colIndices);
    this->mat2Phi( rowIndices, colIndices ).fill( V );
}

void poissonPot::generateGrad2LambdaMatrix( ){
    int i, j, pos;
    double hxMinus1, hxPlus1, hyMinus1, hyPlus1, A, tempEps, EpsNeigh;
    for( j=0; j < this->nNY; j++ ){
        for( i=0; i < this->nNX; i++ ){
            // row or column position in the equation matrix
            pos = i + (j-1) * this->nNX;
            // condition for dirichlet boundary
            if( this->mat2PotentialDirichlet( i, j ) != this->Inf ){
                this->mat2SparseGrad2Lambda( pos, pos ) = 1;
                continue;
            }
            // grid dependent parameters
            this->calculateh(i, j, hxMinus1, hxPlus1, hyMinus1, hyPlus1);
            // other than dirichlet condition
            ////// X Dimension
            if( this->nNX > 1 ){
                // delx i-1
                if( i > 0 ){
                    A = 1 / (  hxMinus1 * ( hxMinus1/2 + hxPlus1/2 )  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i-1, j );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->mat2SparseGrad2Lambda( pos, pos-1 ) = A * tempEps;
                    this->mat2SparseGrad2Lambda(pos, pos) -= A * tempEps;
                }
                // delx i+1
                if( i+1 < this->nNX ){
                    A = 1 / (  hxPlus1 * ( hxMinus1/2 + hxPlus1/2 )  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i+1, j );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->mat2SparseGrad2Lambda( pos, pos+1 ) = A * tempEps;
                    this->mat2SparseGrad2Lambda( pos, pos ) -= A * tempEps;
                }
            }
            ///// Y Dimension
            if( this->nNY > 1 ){
                // dely j-1
                if( j > 0 ){
                    A = 1 / (  hyMinus1 * ( hyMinus1/2 + hyPlus1/2 )  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i, j-1 );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->mat2SparseGrad2Lambda( pos, pos - this->nNX ) = A * tempEps;
                    this->mat2SparseGrad2Lambda( pos, pos ) -= A * tempEps;
                }
                // dely j+1
                if( j+1 < this->nNY ){
                    A = 1 / (  hyPlus1 * ( hyMinus1/2 + hyPlus1/2 )  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i, j+1 );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->mat2SparseGrad2Lambda( pos, pos + this->nNX ) = A * tempEps;
                    this->mat2SparseGrad2Lambda( pos, pos ) -= A * tempEps;
                }
            }
            
        }
    } 
}

mat poissonPot::setInitialGuess( ){
    mat tempmat2Phi;
    vec tempLambda;
    int i, j, pos;
    bool statusSolver;
    for( j=0; j<this->nNY; j++ ){
        for( i=0; i<this->nNX; i++ ){
            pos = i + (j-1) * this->nNX;
            if( this->mat2PotentialDirichlet( i, j ) != this->Inf ){
                this->vecRHS( pos ) = this->mat2Phi( i, j );
            }
        }
    }
    statusSolver = arma::spsolve(tempLambda, this->mat2SparseGrad2Lambda, this->vecRHS);  // use default solver
    //tempLambda = arma::spsolve( this->mat2SparseGrad2Lambda, this->vecRHS);
    
    if( statusSolver == false ){
        cout << "WARNING : Sparse solver found no solution" << endl; 
    }
    tempmat2Phi = mat( tempLambda );
    tempmat2Phi.reshape( this->nNX, this->nNY );
    return tempmat2Phi;
}

mat poissonPot::calculateLambdaSingleIteration( ){
    int i, j, pos;
    double hxMinus1, hxPlus1, hyMinus1, hyPlus1, A, tempEps, EpsNeigh, f1, f2;
    vec tempLambda;
    bool statusSolver;
    mat tempMat2Lambda;
    
    this->vecd2V_by_dx2.zeros();
    this->vecd2V_by_dy2.zeros();
    this->vecRho.zeros();
    this->vecRHS.zeros();
    this->vecdRho_dV.zeros();
    
    for( j=0; j<this->nNX; j++ ){
        for( i=0; i<this->nNY; i++ ){
           pos = i + (j-1) * this->nNX;
           // condition for dirichlet boundary
            if( this->mat2PotentialDirichlet( i, j ) != this->Inf ){
                continue;
            }
            // grid dependent parameters
            this->calculateh(i, j, hxMinus1, hxPlus1, hyMinus1, hyPlus1);
            // X Dimension 
            if( this->nNX > 1 ){
                // delx i-1
                if( i > 0 ){
                    A = 1 / (  hxMinus1 * ( hxMinus1 + hxPlus1 ) / 2  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i-1 , j );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->vecd2V_by_dx2( pos ) += A * tempEps * this->mat2Phi( i-1, j ) - A * tempEps * this->mat2Phi( i, j );
                }
                // delx i+1
                if( i+1 < this->nNX ){
                    A = 1 / (  hxPlus1 * ( hxMinus1 + hxPlus1 ) / 2  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i+1 , j );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->vecd2V_by_dx2( pos ) += A * tempEps * this->mat2Phi( i+1, j ) - A * tempEps * this->mat2Phi( i, j );
                }
            }
            if( this->nNY > 1 ){
                // dely j-1
                if( j > 0 ){
                    A = 1 / (  hyMinus1 * ( hyMinus1 + hyPlus1 ) / 2  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i , j-1 );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->vecd2V_by_dy2( pos ) += A * tempEps * this->mat2Phi( i, j-1 ) - A * tempEps * this->mat2Phi( i, j );
                }
                // dely j+1
                if( j+1 < this->nNY ){
                    A = 1 / (  hyPlus1 * ( hyMinus1 + hyPlus1 ) / 2  );
                    tempEps = this->mat2epsilon( i, j );
                    EpsNeigh = this->mat2epsilon( i , j+1 );
                    if( tempEps != EpsNeigh ){
                        tempEps = tempEps = tempEps * EpsNeigh / (   weightEps * EpsNeigh + ( 1 - weightEps ) * tempEps   );
                    }
                    this->vecd2V_by_dy2( pos ) += A * tempEps * this->mat2Phi( i, j+1 ) - A * tempEps * this->mat2Phi( i, j );
                
                }
                
            }
            // terms for RHS
            if ( (  this->nNX > 1 && ( i > 0 && i+1 < this->nNX)  ) || (  this->nNY > 1 && ( j > 0 && j+1 < this->nNY )   )) {
                double tempV = this->mat2Phi(i, j);
                this->vecRho( pos, 1 ) = this->getRho(i, j, tempV);
                f1 = this->getRho( i, j, tempV + this->DELPHI );
                f2 = this->getRho( i, j, tempV - this->DELPHI );
                this->vecdRho_dV( pos ) = ( f1 - f2 ) / ( 2 * this->DELPHI );
            }
        }
    }
    this->vecRHS = - ( this->vecd2V_by_dx2 + this->vecd2V_by_dy2 + this->vecRho );
    spmat tempMat2SparseGrad2Lambda = this->mat2SparseGrad2Lambda;
    tempMat2SparseGrad2Lambda.diag(0) += this->vecdRho_dV;
    
    statusSolver = arma::spsolve(tempLambda, tempMat2SparseGrad2Lambda, this->vecRHS);  // use default solver
   
    if( statusSolver == false ){
        cout << "WARNING : Sparse solver found no solution" << endl; 
    }
    tempMat2Lambda = mat( tempLambda );
    tempMat2Lambda.reshape( this->nNX, this->nNY );
    return tempMat2Lambda;
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

void poissonPot::calculateh( double xi, double yj, double &hxMinus, double &hxPlus, double &hyMinus, double &hyPlus ){
    
}

void poissonPot::helperGenRowColIndices( double x1, double x2, double y1, double y2, uwcol &rowIndices, uwcol &colIndices){
    uwcol tempXup = find ( this->vecX > x1 );
    uwcol tempXLow = find ( this->vecX <= y1 );
    maths::vintersection( tempXup, tempXLow, rowIndices );
    
    uwcol tempYup = find ( this->vecY > y1 );
    uwcol tempYLow = find ( this->vecY <= y2 );
    maths::vintersection( tempYup, tempYLow, colIndices );
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