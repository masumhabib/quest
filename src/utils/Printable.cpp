/* 
 * File:   Printable.cpp
 * Author: K M Masum Habib <khabib@ee.ucr.edu>
 *
 * Created on April 18, 2013, 4:27 PM
 * 
 */

#include "Printable.hpp"
#include "../python/boostpython.hpp"


/**
 * Python exporter.
 */
namespace qmicad{
namespace python{


void export_Printable(){
    using utils::Printable;
    using utils::stds::string;

    class_<Printable, shared_ptr<Printable> >("Printable", 
            init<optional<const string&> >())
        .def(self_ns::str(self_ns::self))
    ;    

}


}
}



