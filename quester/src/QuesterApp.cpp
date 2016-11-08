/**
 *
 *
 */

#include "QuesterApp.hpp"


namespace quester {

QuesterApp::QuesterApp (int argc, char** argv) {
    //init (argc, argv);
}

void QuesterApp::run() {
}


void QuesterApp::parse_args (int argc, char** argv) {
    //namespace po = boost::program_options;

    //po::options_description generic("General");
    //generic.add_options()
    //    ("version",
    //            "Print version string")
    //    ("help",        
    //            "Produce help message")    
    //    ;
    //    
    //bool batch_mode = opts.mode == Options::Mode::Batch;
    //bool interactive_mode = opts.mode == Options::Mode::Interactive;

    //po::options_description config("Configuration");
    //config.add_options()
    //    ("input-file,i", po::value<std::string>(opts.input_file), 
    //         "Input file name.")
    //    ("batch,B", po::value<bool>(&batch_mode)->default_value(batch_mode), 
    //            "Run in batch mode.")
    //    ("interactive,I", po::value<bool>(&interactive_mode)->default_value(interactive_mode), 
    //            "Run in interactive mode -- much like MATLAB command mode.")
    //    ;
    //
    //// Hidden options, will be allowed both on command line and
    //// in config file, but will not be shown to the user.
    //po::options_description hidden("Hidden options");
    //hidden.add_options()
    //    ("long-version", 
    //            "Print the version info in long, verbose format.")
    //    ;        
}

}

