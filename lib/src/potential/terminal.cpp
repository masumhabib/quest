/* 
 * File:   teminal.cpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
 *
 * Created on January 27, 2014, 4:02 PM
 */

#include "potential/terminal.h"

namespace quest{
namespace potential{

// Construct Terminal from a string containing definition of a plygon
// "POLYGON((x1 y1, x2 y2, ...))"
Terminal::Terminal(const string& polygon, const string &prefix): 
        Printable(" " + prefix)
{
    mTitle = "Terminal";
    bg::read_wkt(polygon, geom);
}

// Construct from a vector of points
Terminal::Terminal(const vector<point>& points, const string &prefix):
        Printable(" " + prefix)
{
    mTitle = "Terminal";
    for (vector<point>::const_iterator it = points.begin(); it != points.end(); ++it){
        geom.outer().push_back(*it);
    }
}

// Construct from a vector of points
Terminal::Terminal(const string &prefix): Printable("" + prefix)
{
    mTitle = "Terminal";
}

string Terminal::toString() const { 
    stringstream ss;
    ss << Printable::toString() << ":" << endl;
    ss << " " << mPrefix << " " << bg::wkt<polygon>(geom);
    return ss.str(); 
}; 

bool Terminal::contains(double x, double y){
    point p(x,y);
    return bg::within(p, geom, stwithin());
}

/*
 * Quadrilateral terminal
 */
Terminal4::Terminal4(point lb, point rb, point rt, point lt, 
        const string& prefix)
{
    Terminal::Prefix(" " + prefix);
    
    geom.outer().push_back(lb);
    geom.outer().push_back(rb);
    geom.outer().push_back(rt);
    geom.outer().push_back(lt);
    
    bg::correct(geom);
}

string Gate4::toString() const{
    stringstream ss;
    ss << Terminal4::toString() << endl;
    ss << " " << mPrefix << " V = " << V;
    
    return ss.str();
}

string Contact4::toString() const{
    stringstream ss;
    ss << Terminal4::toString() << endl;
    ss << " " << mPrefix << " V = " << V << " E = " << E;
    
    return ss.str();
}

}
}
