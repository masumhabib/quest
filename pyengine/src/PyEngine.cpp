/**
 *
 *
 *
 */


#include "PyEngine.hpp"
#include <iostream>


namespace quester {

PyEngine::PyEngine () {
    try {
    	Py_Initialize();
    	main_module = bp::import ("__main__");
    	main_namespace = main_module.attr ("__dict__");
	} catch(bp::error_already_set const&) {
        PyErr_Print();
        //std::string perror_str = bp::parse_python_exception();
        //std::cout << "Error in Python: " << perror_str << std::endl;
    }

}

PyEngine::~PyEngine () {
}

PyEngine::CommandStatus PyEngine::eval (const std::string& command) {
    try {
    	auto result = bp::exec (bp::str (command), main_namespace); 
	} catch(bp::error_already_set const&) {
        PyErr_Print();
    }
    
    return CommandStatus::SUCCESS;
}


void PyEngine::clear () {
}

void PyEngine::clear_history () {
}

std::string PyEngine::get_output_msg () const {
    return std::string ("Hello World!");
}

std::string PyEngine::get_error_msg () const {
    return std::string ();
}


}





