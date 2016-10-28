/**
 *
 *
 *
 */


#include "PyEngine.hpp"


namespace quester {

PyEngine::PyEngine () {
}

PyEngine::~PyEngine () {
}

PyEngine::CommandStatus PyEngine::run_command (const std::string& command, 
std::string& output) {
    return CommandStatus::SUCCESS;
}

std::string get_error_msg () {
    return std::string ();
}

void clear () {
}

}





