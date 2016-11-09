/**
 *
 *
 */

#include "QuesterApp.hpp"


namespace quester {

QuesterApp::QuesterApp (int argc, char** argv) {
    parse_args (argc, argv);
}

void QuesterApp::run() {
    if (opts.print_help_message == true) {
        print_usage ();
    }
}


void QuesterApp::parse_args (int argc, char** argv) {
    namespace po = boost::program_options;

    po::options_description general("General");
    general.add_options()
        ("help,h",        
                "Produce help message")
        ("version",
                "Print version string")
        ;
        
    bool batch_mode = false;
    bool interactive_mode = false;
    std::string script_name;

    po::options_description config("Configuration");
    config.add_options()
        ("run-command,r", po::value<std::string>(&opts.command), 
                "Command to run ... e.g., run filename.py")
        ("script_file,f", po::value<std::string>(&script_name), 
                "Script file to run.")
        ("batch,b", po::value<bool>(&batch_mode)->default_value(batch_mode), 
                "Run in batch mode.")
        ("interactive,i", po::value<bool>(&interactive_mode)->default_value(interactive_mode), 
                "Run in interactive mode -- much like MATLAB command mode.")
        ;

    // Hidden options, will be allowed both on command line and
    // in config file, but will not be shown to the user.
    po::options_description hidden("Hidden options");
    hidden.add_options()
        ("long-version", 
                "Print the version info in long, verbose format.")
        ;        

    // command line options
    po::options_description cmd_line_opts;
    cmd_line_opts.add (general).add (config).add (hidden);

    // config file options
    po::options_description conf_file_opts;
    conf_file_opts.add (general).add (config).add (hidden);

    // shown in the help message
    po::options_description usage ("");
    usage.add (general).add (config);
    std::stringstream out;
    out << usage;
    usage_message = out.str();

    // do the actual parsing
    try {
        po::variables_map users_choices;
        po::store (po::command_line_parser (argc, argv)
                        .options (cmd_line_opts).run(), users_choices);
        po::notify (users_choices);


        //if (!batch_mode && !interactive_mode) {
        //    if (opts.script_name.length() == 0) {
        //        opts.mode = Options::Mode::Interactive;
        //    }
        //}
        

        if (users_choices.count ("help")) {
            opts.print_help_message = true;
        }

    } catch (const po::unknown_option & e) {
        usage_error_message = e.what();
        opts.print_help_message = true;
    }
}

void QuesterApp::print_usage () {
    if (usage_error_message.length () > 0) {
        std::cout << std::endl << "ERROR: " << usage_error_message 
                  << std::endl << std::endl;
    }
    std::cout << usage_message << std::endl;
}

}

