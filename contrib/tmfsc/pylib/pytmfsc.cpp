/***/
#include "pydevice.h"
#include "boostpython.hpp"
#include "qmicad.hpp"



namespace qmicad{
namespace python{

using namespace boost::python;


BOOST_PYTHON_MODULE(tmfsc)
{ 
    // create qmicad package
    //object package = scope();
    //package.attr("__path__") = "tmfsc";
    //package.attr("version") = version;

    //def("greet", greet, " Shows the QMICAD banner.");
    //def("setVerbosity", setVerbosity, " Sets the verbosity level of C++ code.");

    export_Device();
    
}

}}

