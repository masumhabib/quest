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

    string (Printable::*Printable_Title1)() = &Printable::Title;
    void (Printable::*Printable_Title2)(const string &title) = &Printable::Title;
    string (Printable::*Printable_Prefix1)() = &Printable::Prefix;
    void (Printable::*Printable_Prefix2)(const string& prefix) = &Printable::Prefix;
    class_<Printable, shared_ptr<Printable> >("Printable", 
            init<optional<const string&> >())
        .add_property("Title", Printable_Title1, Printable_Title2)
        .add_property("Prefix", Printable_Prefix1, Printable_Prefix2)
        .def(self_ns::str(self_ns::self))
    ;    

}


}
}



