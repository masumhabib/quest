/* 
 * File:   Lattice.cpp
 * Copyright (C) 2013-2014  K M Masum Habib <masum.habib@gmail.com>
 * 
 * Created on April 7, 2013, 9:46 AM
 * 
 * Description: Lattice vector and lattice coordinates.
 * 
 */

#include "atoms/Lattice.h"
#include "python/boostpython.hpp"

namespace qmicad{
namespace atoms{


LatticeCoordinate::LatticeCoordinate(int n1, int n2, int n3, const string &prefix)
    :Printable(" " + prefix), n1(n1), n2(n2), n3(n3)
{
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

string LatticeCoordinate::toString() const{
    stringstream out;
    
    out << "[";
    out.width(3);
    out << n1; 
    out << ", ";
    out.width(3);
    out << n2; 
    out << ", ";
    out.width(3);
    out << n3;
    out << "]";

    return out.str();
}

/** 
 * Lattice vector.
 */
LatticeVector::LatticeVector(const string &prefix):
    Printable(" " + prefix), a1(3), a2(3), a3(3)
{
    zeros();
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


svec operator*(const LatticeVector &lv, const LatticeCoordinate& lc){
    svec sv = lv.a1*lc.n1 + lv.a2*lc.n2 + lv.a3*lc.n3;
    return sv;
}

string LatticeVector::toString() const{
        
    stringstream out;
    out.precision(4);
    out << " a1 = [" << std::fixed << a1(coord::X) << " " 
                     << std::fixed << a1(coord::Y) << " " 
                     << std::fixed << a1(coord::Z) << "]" << endl;
    out << " a2 = [" << std::fixed << a2(coord::X) << " " 
                     << std::fixed << a2(coord::Y) << " " 
                     << std::fixed << a2(coord::Z) << "]" << endl;
    out << " a3 = [" << std::fixed << a3(coord::X) << " " 
                     << std::fixed << a3(coord::Y) << " " 
                     << std::fixed << a3(coord::Z) << "]";
    
    return out.str();
}

}
}


/**
 * Python exporter.
 */
namespace qmicad{
namespace python{

void export_lvec(){
    using namespace atoms;
    
    /**
     * Lattice vector.
     */
    class_<lvec, bases<Printable>, shared_ptr<lvec> >("LVec", 
            init<optional<const string&> >())
        .def_readwrite("a1", &lvec::a1)
        .def_readwrite("a2", &lvec::a2)
        .def_readwrite("a3", &lvec::a3)
        .def("la1", &lvec::la1)
        .def("la2", &lvec::la2)
        .def("la3", &lvec::la3)
    ;
}


void export_lcoord(){
    using namespace atoms;
    /**
     * Lattice coordinate.
     */
    class_<lcoord, bases<Printable>, shared_ptr<lcoord> >("LCoord", 
            init<int, int, int, optional<const string&> >())
        .def_readwrite("n1", &lcoord::n1)
        .def_readwrite("n2", &lcoord::n2)
        .def_readwrite("n3", &lcoord::n3)
    ;
}


}
}


