
/* 
 * File:   read_python.h
 * Copyright (C) 2014 K M Masum Habib <masum.habib@ail.com>
 *
 * Created on November 3, 2014, 2:59 PM
 */

#ifndef READ_PYTHON_H
#define	READ_PYTHON_H

#include <boostpython.hpp>
#include <boost/python.hpp>
#include <boost/python/stl_iterator.hpp>
#include <string>
#include <sstream>
#include <vector>

namespace utils{ namespace python{

    namespace bp = boost::python;
    using std::string;
    using std::vector;
    using std::stringstream;

    void create_main_module(bp::object &main_module);
    void eval_options(const string& file_name, bp::object &module);
    bool hasattr(const bp::object& module, const string& name);
    
    /**
     * Reads a scalar option from python script.
     */
    template <typename T>
    bool read_option(const bp::object& module, const string& name, T& value){
        if (hasattr(module, name)){
            bp::object dict = module.attr("__dict__");
            value = bp::extract<T>(dict[name]);
            return true;
        }    
        return false;
    }

    /**
     * Reads a vector option from python script.
     */
    template<typename T>
    bool read_option(const bp::object& module, const string& name, vector<T>& value) {
        if (hasattr(module, name)){
            bp::object dict = module.attr("__dict__");
            bp::stl_input_iterator<T> begin(dict[name]);
            bp::stl_input_iterator<T> end;
            value.clear();
            value.insert(value.end(), begin, end); 
            return true;
        } 
        return false;
    }

    /**
     * Writes a scalar to a python script.
     */
    template <typename T>
    bool write_option(bp::object module, const string& name, const T& value){
        stringstream py_cmd;
        py_cmd << name << " = " << value;
        bp::object nmspace = module.attr("__dict__");
        bp::object ignored = bp::exec(bp::str(py_cmd.str()), nmspace); 
        return true;
    }

    
    
}}
#endif	/* READ_PYTHON_H */

