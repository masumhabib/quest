/**
 *
 *
 *
 */

#ifndef QUESTERAPP_HPP
#define QUESTERAPP_HPP

#include "qlog.hpp"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <sstream>
#include <fstream>

namespace quester {

struct Options {
public:
    enum class Mode {
        None = 0, Batch, Interactive,
    };


    Mode mode = Mode::None;
    std::string command;

    bool print_help_message = false;

    qlog::Verbosity verbosity;
};

class QuesterApp {
public:
    QuesterApp (int argc, char** argv);
    void run();

private:
    void init ();
    void parse_args (int argc, char** argv);
    void print_usage ();


private:
    Options opts;
    std::string usage_error_message;
    std::string usage_message;
    std::string exe_name;
};

}

#endif



