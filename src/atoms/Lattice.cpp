/* 
 * File:   Lattice.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 * 
 * Created on April 7, 2013, 9:46 AM
 * 
 * Description: Lattice vector and lattice coordinates.
 * 
 */

#include "Lattice.h"

namespace qmicad{

/* Lattice coordinate */

LatticeCoordinate::LatticeCoordinate(){
}

LatticeCoordinate::LatticeCoordinate(int n1, int n2, int n3){
        this->n1 = n1;
        this->n2 = n2;
        this->n3 = n3;
}

LatticeCoordinate::~LatticeCoordinate(){
}

LatticeCoordinate operator+ (LatticeCoordinate lhs, const LatticeCoordinate& rhs){
    lhs += rhs;
    return lhs;
}

LatticeCoordinate& LatticeCoordinate::operator+= (const LatticeCoordinate& rhs){
    
    n1 += rhs.n1;
    n2 += rhs.n2;
    n3 += rhs.n3;
    
    return *this;
}

LatticeCoordinate operator- (LatticeCoordinate lhs, const LatticeCoordinate& rhs){
    lhs -= rhs;
    return lhs;
}

LatticeCoordinate& LatticeCoordinate::operator-= (const LatticeCoordinate& rhs){
    
    n1 -= rhs.n1;
    n2 -= rhs.n2;
    n3 -= rhs.n3;
    
    return *this;
}

ostream& operator << (ostream& out, const LatticeCoordinate& lc){
    
    out << "[";
    out.width(3);
    out << lc.n1; 
    out << ", ";
    out.width(3);
    out << lc.n2; 
    out << ", ";
    out.width(3);
    out << lc.n3;
    out << "]";

}

/* Lattice vector */
LatticeVector::LatticeVector() {
    zeros();
}

LatticeVector::~LatticeVector() {
}

void LatticeVector::zeros(){
    a1.zeros();
    a2.zeros();
    a3.zeros();    
}

LatticeVector operator+ (LatticeVector lhs, const LatticeVector& rhs){
    lhs += rhs;
    return lhs;
}

LatticeVector& LatticeVector::operator+= (const LatticeVector& rhs){
    
    a1 += rhs.a1;
    a2 += rhs.a2;
    a3 += rhs.a3;
    
    return *this;
}

LatticeVector operator- (LatticeVector lhs, const LatticeVector& rhs){
    lhs -= rhs;
    return lhs;
}

LatticeVector& LatticeVector::operator-= (const LatticeVector& rhs){
    
    a1 -= rhs.a1;
    a2 -= rhs.a2;
    a3 -= rhs.a3;
    
    return *this;
}

ostream& operator << (ostream& out, const LatticeVector& lv){
        
    out.precision(4);
    out << " a1 = [" << std::fixed << lv.a1(X) << " " 
                     << std::fixed << lv.a1(Y) << " " 
                     << std::fixed << lv.a1(Z) << "]" << endl;
    out << " a2 = [" << std::fixed << lv.a2(X) << " " 
                     << std::fixed << lv.a2(Y) << " " 
                     << std::fixed << lv.a2(Z) << "]" << endl;
    out << " a3 = [" << std::fixed << lv.a3(X) << " " 
                     << std::fixed << lv.a3(Y) << " " 
                     << std::fixed << lv.a3(Z) << "]";
    
    return out;

}

}

