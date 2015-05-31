/** 
 * Python wrapper for the Device class 
 * author: K M Masum Habib
 */

#include "pydevice.h"

namespace qmicad{
namespace python{

void export_Device(){
    using namespace qmicad::tmfsc;
    class_<Device, bases<Printable>, shared_ptr<Device> >("Device", 
            init<>())
        .def("addPoint", &Device::addPoint)
        .def("addPoints", &Device::addPoints)
    ;
}

}}



