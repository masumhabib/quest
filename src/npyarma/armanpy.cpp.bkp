/**
 *       Filename:  armanpy.cpp
 *
 *    Description:  Converter between armadillo and numpy
 *
 *        Created:  06/23/2014 03:50:57 PM
 *         Author:  K M Masum Habib<masum.habib@gmail.com>
 *      Copyright:  (c) 2014 K M Masum Habib<masum.habib@gmail.com>.
 *
 */

#include "armanpy/armanpy.h"

#include "utils/vout.h"

namespace qmicad{
namespace armanpy{


object vec_to_numpy(){
    using namespace utils::stds;

    vec v(3);
    v << 1 << endr << 2 << endr << 3<<endr;
    cout << v << endl;
    
    //npy_intp size = v.n_cols;
    //PyObject *pyObj = PyArray_SimpleNewFromData(1, &size, NPY_DOUBLE, v.memptr());
    //handle<> handle(pyObj);
    //numeric::array arr(handle);

    //return arr.copy();
}

}
}

namespace qmicad{
namespace python{
using namespace utils::stds;
using namespace armanpy;

void export_armanpy(){
    def("vec_to_numpy", vec_to_numpy, " vec to numpy.array");
}

}
}
