/***/
#include "pytmfsc.h"
#include "quest.hpp"



namespace quest{
namespace python{

using namespace boost::python;


BOOST_PYTHON_MODULE(tmfsc)
{ 
    using namespace quest::tmfsc;
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

