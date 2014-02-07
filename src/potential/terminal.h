/* 
 * File:   teminal.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 27, 2014, 4:02 PM
 * 
 * The terminal class stores terminal properties such as 
 * drain, source and gates
 * 
 */

#ifndef TEMINAL_H
#define	TEMINAL_H

#include <vector>
#include "../utils/Printable.hpp"
#include "../utils/geometry.hpp"

using std::vector;

/*
 * Terminal class stores terminal data structure.
 */
struct Terminal:public Printable{
    polygon geom;       // geometry of the terminal

    // construct from a string containing definition of a plygon
    // "POLYGON((x1 y1, x2 y2, ...))"
    Terminal(const string& polygon, const string &prefix = "");
    // Construct from a vector of points
    Terminal(const vector<point>& points, const string &prefix = "");
    // Without the geometry
    Terminal(const string &prefix = "");
    // Converts to string for << operator
    string toString() const;
    // checks if x,y is within geom
    virtual bool contains(double x, double y);
};

/*
 * Quadrilateral terminal
 */
struct Terminal4:public Terminal{
    // Construct from a four points
    Terminal4(point lb, point rb, point rt, point lt, const string &prefix = "");
};

/*
 * Gates have dirichlet boundaty.
 */
struct Gate4:public Terminal4{
    double  V;         // gate voltage. Dirichlet boundary.
    Gate4(point lb = point(0,0), point rb = point(0,0), 
        point rt = point(0,0), point lt = point(0,0), double V = 0.0, 
        const string &prefix = ""):
        Terminal4(lb, rb, rt, lt, prefix), V(V)
    {
        mTitle = "Quadrilateral Gate";
    }
    virtual string toString() const;
};

/*
 * Contacts have neumann boundary.
 */
struct Contact4:public Terminal4{
    double  E;         // contact electric field. Neumann boundary.
    double  V;         // contact potential. Dirichlet boundary
    Contact4(point lb = point(0,0), point rb = point(0,0), 
        point rt = point(0,0), point lt = point(0,0), double V = 0.0, 
        double E = 0.0, const string &prefix = ""):
        Terminal4(lb, rb, rt, lt, prefix), E(E), V(V)
    {
        mTitle = "Quadrilateral Contact";
    }
    virtual string toString() const;
};

typedef Terminal4   terminal;
typedef Contact4    contact;
typedef Gate4       gate;

#endif	/* TEMINAL_H */

