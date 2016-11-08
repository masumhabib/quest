/**
 *
 *
 *
 */

#ifndef QUESTERAPP_HPP
#define QUESTERAPP_HPP

#include "boost/program_options.hpp" 
//#include "logger/logger.hpp"

namespace quester {

struct Options {
public:
    //enum class Mode {
    //    Batch = 0, Interactive;
    //};


    //Mode mode = Mode::Interactive;
    //std::string input_file;

    //logger::Verbosity verbosity;
};

class QuesterApp {
public:
    QuesterApp (int argc, char** argv);
    void run();

private:
    void init ();
    void parse_args (int argc, char** argv);


private:
    Options opts;
};

}

#endif



