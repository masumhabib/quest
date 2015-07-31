/** 
 * Python wrapper for the Device class 
 * author: K M Masum Habib
 */

#include "pydevice.h"

namespace qmicad{
namespace python{

point PyDevice::contMidPointPy(int indx) {
    int iEdge = contToEdgeIndx(indx);
    return edgeMidPointPy(iEdge);
}

point PyDevice::edgeMidPointPy(int indx) {
    Edge e = edge(indx);
    point midp = {(e.p()[0]+e.q()[0])/2, (e.p()[1]+e.q()[1])/2};
    return midp;
}

void export_Device(){
    using namespace qmicad::tmfsc;

    int (PyDevice::*PyDevice_edgeType1)(int iEdge) = &PyDevice::edgeType;
    void (PyDevice::*PyDevice_edgeType2)(int iEdge, int type) = &PyDevice::edgeType;
    class_<PyDevice, bases<Printable>, shared_ptr<PyDevice> >("Device", 
            init<>())
        .def("addPoint", &PyDevice::addPoint)
        .def("addPoints", &PyDevice::addPoints)
        .def("addEdge", &PyDevice::addEdge)
        .def("addGate", &PyDevice::addGate)
        .def("edgeType", PyDevice_edgeType1)
        .def("edgeType", PyDevice_edgeType2)
        .def("edgeUnitVect", &PyDevice::edgeUnitVectPy)
        .def("edgeNormVect", &PyDevice::edgeNormVect)
        .def("edgeVect", &PyDevice::edgeVect)
        .def("edgeMidPoint", &PyDevice::edgeMidPointPy)
        .def("contMidPoint", &PyDevice::contMidPointPy)
        .def("contNormVect", &PyDevice::contNormVect)
        .def("contDirctn", &PyDevice::contDirctn)
        .def("numEdges", &PyDevice::numEdges)
        .def("numConts", &PyDevice::numConts)
        .def("isReflectEdge", &PyDevice::isReflectEdge)
        .def("isAbsorbEdge", &PyDevice::isAbsorbEdge)
        .def_readonly("NumGates", &PyDevice::getNumGates)
    	.def("getSplitLen", &PyDevice::getSplitLen)
    	.def("setSplitLen", &PyDevice::setSplitLen);
 
    ;
}

}}



