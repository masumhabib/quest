/** 
 * Python wrapper for the Device class 
 * author: K M Masum Habib
 */

#include "pydevice.h"

namespace qmicad{
namespace python{

void export_Device(){
    using namespace qmicad::tmfsc;

    int (PyDevice::*PyDevice_edgeType1)(int iEdge) = &PyDevice::edgeType;
    void (PyDevice::*PyDevice_edgeType2)(int iEdge, int type) = &PyDevice::edgeType;
    class_<PyDevice, bases<Printable>, shared_ptr<PyDevice> >("Device", 
            init<>())
        .def("addPoint", &PyDevice::addPoint)
        .def("addPoints", &PyDevice::addPoints)
        .def("edgeType", PyDevice_edgeType1)
        .def("edgeType", PyDevice_edgeType2)
        .def("numEdges", &PyDevice::numEdges)
        .def("isReflectEdge", &PyDevice::isReflectEdge)
        .def("isAbsorbEdge", &PyDevice::isAbsorbEdge)
        .def("edgeUnitVect", &PyDevice::edgeUnitVectPy)
        .def("edgeVect", &PyDevice::edgeVect)
    ;
}

}}



