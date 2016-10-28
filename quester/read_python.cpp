
/* 
 * File:   read_python.cpp
 * Copyright (C) 2014 K M Masum Habib <masum.habib@ail.com>
 *
 * Created on November 3, 2014, 2:59 PM
 */

#include "read_python.h"

namespace utils{ namespace python{

void create_main_module(bp::object &main_module){
    main_module = bp::import("__main__");
}

void eval_options(const string& file_name, bp::object &module){
    bp::object nmspace = module.attr("__dict__");
    bp::object ignored = bp::exec_file(bp::str(file_name), nmspace);
}

// Helper function to check if an object has an attribute.
bool hasattr(const bp::object& module, const string& name){
    return PyObject_HasAttrString(module.ptr(), name.c_str());
}


}}



