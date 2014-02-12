/* 
 * File:   PyVecGrid.h
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on February 11, 2014, 1:08 AM
 */

#ifndef PYVECGRID_H
#define	PYVECGRID_H

#include "../grid/grid.hpp"

namespace qmicad{
namespace python{
using utils::VecGrid;
using utils::stds::string;

class PyVecGrid: public VecGrid{
public:
    PyVecGrid(double min, double max, double d, const string &prefix = "")
            :VecGrid(min, max, d, prefix)
    {            
    }
    
    PyVecGrid(double min = 0, double max = 0, int N = 1, const string &prefix = ""):
            VecGrid(min, max, N, prefix)
    {    
    }
    
    double V(int it){
        return (*this)(it);
    }
};

}
}
#endif	/* PYVECGRID_H */

