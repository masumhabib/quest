/**
 * @brief This class embeds python engine and provides high level interface 
 *        to C++.
 *
 *
 */


#ifndef QUESTER_PYENGINE_HPP
#define QUESTER_PYENGINE_HPP

#include <string>
#include <vector>
#include <boost/python.hpp>


namespace quester {
namespace bp = boost::python;

class PyEngine {
public:
    enum class CommandStatus {
        ERROR   = 0,
        SUCCESS = 1,
    };

    struct Options {
    };

public:
    PyEngine ();
    PyEngine (const Options& opts);
    ~PyEngine ();

    CommandStatus eval (const std::string& command);

    void clear ();
    void clear_history ();

    std::vector <std::string> get_history () const;
    std::string get_output_msg () const;
    std::string get_error_msg () const;

private:
    bp::object main_module;
    bp::object main_namespace;
};

}

#endif








