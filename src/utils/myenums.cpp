/* 
 * File:   myenums.hpp
 * Copyright (C) 2014  K M Masum Habib <masum.habib@gmail.com>
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



