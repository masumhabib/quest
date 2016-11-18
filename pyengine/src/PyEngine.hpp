/**
 *
 *
 *
 */


#ifndef QUESTER_PYENGINE_HPP
#define QUESTER_PYENGINE_HPP

#include <string>


namespace quester {

class PyEngine {
public:
    enum class CommandStatus {
        ERROR   = 0,
        SUCCESS = 1,
    };

public:
    PyEngine ();
    ~PyEngine ();

    CommandStatus run_command (const std::string& command, std::string& output);

    std::string get_error_msg ();
    void clear ();
};

}

#endif








