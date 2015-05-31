/***/
#include "pytmfsc.h"
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

    export_Device();
    export_Simulator();
    
}

}}

