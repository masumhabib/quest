/* 
 * File:   myenums.hpp
 * Author: K M Masum Habib<masum.habib@virginia.edu>
 *
 * Created on January 21, 2014, 11:42 AM
 */

#include "myenums.hpp"
#include "../python/boostpython.hpp"

namespace qmicad{
namespace python{

/**
* python export.
*/
void export_Option(){
    using namespace utils::enums;
    
    enum_<Option>("Option") 
       .value("Disabled", Disabled)
       .value("Enabled",  Enabled)
    ;
}

}
}



