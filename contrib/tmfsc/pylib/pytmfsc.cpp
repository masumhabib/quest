/***/
#include "pytmfsc.h"
#include "boostpython.hpp"
#include "qmicad.hpp"



namespace qmicad{
namespace python{

using namespace boost::python;


BOOST_PYTHON_MODULE(tmfsc)
{ 
    using namespace qmicad::tmfsc;
    // add global variables
    object package = scope();
    package.attr("nm") = nm;
    package.attr("AA") = AA;
    package.attr("EDGE_ABSORB") = Edge::EDGE_ABSORB;
    package.attr("EDGE_REFLECT") = Edge::EDGE_REFLECT;
    package.attr("EDGE_TRANSMIT") = Edge::EDGE_TRANSMIT;

    export_Device();
    export_Simulator();
    
}

}}

